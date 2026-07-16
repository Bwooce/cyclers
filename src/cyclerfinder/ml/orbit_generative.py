"""#608 proof-of-concept: statistical generative seed model for CR3BP periodic orbits.

Externally de-risked by Litteri, Gil, Vasile, Rodriguez-Fernandez & Camacho,
"Generation of periodic orbits in the restricted three-body problem with a
variational autoencoder," *Celestial Mechanics and Dynamical Astronomy* 138:25
(June 2026) -- a CNN-VAE trained on ~44k time-series-represented CR3BP orbits
across 40 Earth-Moon families, sampled from the 2D latent space, and refined
via multiple-shooting (46% of 100 latent samples converged to a genuinely new
orbit). See that paper's earlier conference version (arXiv:2408.03691) for the
architecture details this module was checked against.

**Design decision (bounded POC, no new ML dependency).** This project has no
torch/jax dependency (see ``pyproject.toml``). A literal CNN-VAE is hard to
implement well with numpy/scipy alone, and this task is explicitly scoped as a
bounded feasibility check, not a production pipeline -- adding a deep-learning
framework for a single proof-of-concept run does not clear this project's
"don't add abstractions beyond what's needed" bar. Litteri et al.'s CNN
encoder/decoder operates on a 100-node time-series representation of each
orbit ONLY because their architecture needs a fixed-length sequence to exploit
temporal structure; the physics does not require it. A CR3BP periodic orbit is
already exactly and minimally described by its 6-component initial state plus
period (an 8-vector once the derived Jacobi constant is appended for
convenience) -- that IS the "reduced coordinate space" the task's own framing
invites as a lighter alternative. ``ClusteredGaussianLatentModel`` below is
the linear-Gaussian analog of their approach: standardize -> PCA (linear
encoder/decoder, replacing the CNN) -> k-means partition -> per-cluster
empirical Gaussian (replacing the learned, single isotropic-Gaussian latent
prior with a piecewise one, since #608's assembled corpus is visibly
multi-modal -- distinct continuation branches/families, not one blob). It
tests the SAME hypothesis (can a statistical model of "what a valid converged
orbit's genome looks like" generate useful NEW seeds) without the heavier
dependency. If this bounded test shows real promise, revisit whether a true
VAE is warranted -- do not build one speculatively first.

Training data provenance: this module reads the project's OWN accumulated
``#210`` outcome-log corpus (``out/outcome_log/*.jsonl``, produced by
``cyclerfinder.search.outcome_log.log_outcome`` as a passive byproduct of past
corrector runs -- never synthetic/toy data). That corpus is large (~540k raw
lines from a dozen past campaigns) but heterogeneous: many "converged" records
are numerically converged to degenerate, non-physical solutions (near-collision
tiny-period loops, huge-velocity junk swept up by continuation excursions in
mu/x0). :func:`assemble_corpus` applies explicit physical-plausibility filters
(Jacobi/period/out-of-plane-amplitude bounds) documented inline below; see
``scripts/run_608_generative_seed_poc.py`` and its report for the histogram
that motivated the specific cutoffs.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.cluster.vq import kmeans2
from scipy.optimize import minimize
from scipy.spatial import cKDTree

FEATURE_NAMES: tuple[str, ...] = ("x0", "y0", "z0", "vx0", "vy0", "vz0", "period", "jacobi")

# Default physical-plausibility window for an Earth-Moon CR3BP periodic orbit
# (see the module docstring: the raw #210 outcome log's "converged" flag alone
# lets through many degenerate solves -- near-collision tiny-period loops,
# huge-velocity junk from continuation excursions -- that are numerically
# converged but not physically meaningful periodic orbits). Shared by
# ``assemble_corpus`` (filters the TRAINING corpus) and ``is_physically_sane``
# (classifies a CORRECTOR OUTPUT after refinement) so both sides of the #608
# generate-then-refine comparison use the identical yardstick.
DEFAULT_JACOBI_BOUNDS: tuple[float, float] = (2.5, 3.5)
DEFAULT_PERIOD_BOUNDS: tuple[float, float] = (0.5, 15.0)
DEFAULT_MAX_ABS_Z0: float = 0.3
DEFAULT_MAX_ABS_VZ0: float = 1.0


def is_physically_sane(
    state0: NDArray[np.float64] | list[float],
    period: float,
    jacobi: float,
    *,
    jacobi_bounds: tuple[float, float] = DEFAULT_JACOBI_BOUNDS,
    period_bounds: tuple[float, float] = DEFAULT_PERIOD_BOUNDS,
    max_abs_z0: float = DEFAULT_MAX_ABS_Z0,
    max_abs_vz0: float = DEFAULT_MAX_ABS_VZ0,
) -> bool:
    """Classify a (possibly corrector-refined) orbit as physically plausible.

    A solver's ``converged=True`` flag means only "the Newton iteration drove
    the periodicity residual below tolerance" -- it says nothing about
    whether the result is a physically meaningful member of a real Earth-Moon
    CR3BP family rather than a degenerate/collision-adjacent curiosity. This
    is the SAME bound used to build the training corpus (see
    ``assemble_corpus``), applied post-hoc to a corrector's output so
    "did refinement succeed" can be reported honestly as "succeeded AND
    landed somewhere physically sane," not just "the residual was small."
    """
    if not (jacobi_bounds[0] <= jacobi <= jacobi_bounds[1]):
        return False
    if not (period_bounds[0] <= period <= period_bounds[1]):
        return False
    return abs(state0[2]) <= max_abs_z0 and abs(state0[5]) <= max_abs_vz0


def heuristic_family_tag(
    state0: NDArray[np.float64] | list[float],
    period: float,
    jacobi: float,
    *,
    jacobi_bounds: tuple[float, float] = DEFAULT_JACOBI_BOUNDS,
    period_bounds: tuple[float, float] = DEFAULT_PERIOD_BOUNDS,
    n_jacobi_bins: int = 3,
    n_period_bins: int = 3,
    planar_z_threshold: float = 1e-3,
) -> str:
    """Cheap, HONEST heuristic family-*proxy* tag -- #614's answer to `#608`'s

    flagged "no family label" gap. This is explicitly NOT a validated
    halo/Lyapunov/DRO/NRHO orbit-family classification: no such classifier
    for an arbitrary ``(state0, period, jacobi)`` genome exists anywhere in
    this codebase (checked ``src/cyclerfinder/search/`` and
    ``src/cyclerfinder/genome/`` -- the "family" hits there are either
    sourced literature *anchors* for specific known families, e.g.
    ``halo_family_at_jacobi.py``/``genome/known_corpus_3d.py``, or family-
    SWITCHING/continuation machinery, e.g. ``genome/family_switch.py``, none
    of which classifies an arbitrary already-converged orbit's family).

    Instead, bins each orbit by three quantities that correlate strongly
    with which real CR3BP family/branch it plausibly belongs to -- energy
    (Jacobi band) roughly sets which families exist at a given level, period
    (log-spaced band, since this corpus spans an order of magnitude)
    distinguishes branches/resonances within a family, and out-of-plane
    sign/magnitude (``z0`` sign, with a small "planar" deadband) separates
    northern/southern halo-like branches from planar (Lyapunov/DRO-like)
    orbits. This does NOT amount to a rigorous classification -- two
    genuinely different real families can land in the same
    Jacobi/period/z-sign bin and get the same tag, and one real family can
    straddle a bin boundary and get split across two tags. Present this tag
    honestly as a heuristic proxy, never as a validated family label.
    """
    j_lo, j_hi = jacobi_bounds
    p_lo, p_hi = period_bounds
    j_clamped = min(max(float(jacobi), j_lo), j_hi)
    p_clamped = min(max(float(period), p_lo), p_hi)
    j_bin = int(
        np.clip(np.floor((j_clamped - j_lo) / (j_hi - j_lo) * n_jacobi_bins), 0, n_jacobi_bins - 1)
    )
    # Period spans roughly an order of magnitude in this corpus -> log-spaced bins.
    log_p_lo, log_p_hi = np.log(p_lo), np.log(p_hi)
    log_p = np.log(max(p_clamped, 1e-9))
    p_bin = int(
        np.clip(
            np.floor((log_p - log_p_lo) / (log_p_hi - log_p_lo) * n_period_bins),
            0,
            n_period_bins - 1,
        )
    )
    z0 = float(state0[2])
    if abs(z0) < planar_z_threshold:
        z_tag = "z0"
    elif z0 > 0:
        z_tag = "z+"
    else:
        z_tag = "z-"
    return f"J{j_bin}_P{p_bin}_{z_tag}"


@dataclass(frozen=True)
class OrbitCorpus:
    """Assembled, deduplicated ``(state0, period, jacobi)`` records for one CR3BP pair."""

    features: NDArray[np.float64]  # (n, 8): x0,y0,z0,vx0,vy0,vz0,period,jacobi
    mu: float
    primary: str
    secondary: str
    n_scanned: int
    n_converged_prefilter: int
    # #614: cheap heuristic family-proxy tag per row (see `heuristic_family_tag`
    # docstring for exactly what this is/isn't). ``None`` only for corpora built
    # before this field existed; `assemble_corpus` always populates it.
    family_tags: NDArray[np.str_] | None = None

    def __len__(self) -> int:
        return int(self.features.shape[0])


def iter_outcome_records(paths: Iterable[Path]) -> Iterator[dict]:
    """Yield parsed JSON objects from #210 outcome-log JSONL files.

    Silently skips lines that fail to parse (the log is append-only,
    best-effort, and occasionally interrupted mid-write by a killed worker
    process -- a torn last line must never abort assembly of the other
    hundreds of thousands of good ones).
    """
    for path in paths:
        with Path(path).open(encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue


def assemble_corpus(
    paths: Iterable[Path],
    *,
    solver: str = "cr3bp.correct_periodic",
    primary: str = "Earth",
    secondary: str = "Moon",
    jacobi_bounds: tuple[float, float] = DEFAULT_JACOBI_BOUNDS,
    period_bounds: tuple[float, float] = DEFAULT_PERIOD_BOUNDS,
    max_residual: float = 1e-9,
    max_abs_z0: float = DEFAULT_MAX_ABS_Z0,
    max_abs_vz0: float = DEFAULT_MAX_ABS_VZ0,
    dedup_decimals: int = 4,
) -> OrbitCorpus:
    """Assemble a REAL (not synthetic) training corpus from #210 outcome logs.

    Keeps only records that are genuinely converged ``correct_periodic``
    solves for ``(primary, secondary)`` AND fall within a physically
    plausible window. The raw log's "converged" flag alone is not a
    sufficient physical-plausibility filter -- see this module's docstring.
    Near-duplicate rows (dense continuation stepping revisits nearly the same
    state many times) are collapsed by rounding to ``dedup_decimals`` places
    before set-based dedup, so the assembled corpus reflects distinct genomes
    rather than an inflated count along one curve.

    Raises
    ------
    ValueError
        If no record survives the filters (an honest failure, not a silent
        empty corpus).
    """
    n_scanned = 0
    n_converged = 0
    seen: set[tuple[float, ...]] = set()
    rows: list[tuple[float, ...]] = []
    tags: list[str] = []
    for rec in iter_outcome_records(paths):
        n_scanned += 1
        if rec.get("solver") != solver:
            continue
        meta = rec.get("meta") or {}
        if meta.get("primary") != primary or meta.get("secondary") != secondary:
            continue
        outcome = rec.get("outcome") or {}
        if not outcome.get("converged"):
            continue
        n_converged += 1
        residual = outcome.get("residual")
        jacobi = outcome.get("jacobi")
        period = outcome.get("period")
        if residual is None or jacobi is None or period is None:
            continue
        if residual >= max_residual:
            continue
        if not (jacobi_bounds[0] <= jacobi <= jacobi_bounds[1]):
            continue
        if not (period_bounds[0] <= period <= period_bounds[1]):
            continue
        state0 = rec.get("inputs", {}).get("state0_guess")
        if not state0 or len(state0) != 6:
            continue
        if abs(state0[2]) > max_abs_z0 or abs(state0[5]) > max_abs_vz0:
            continue
        row = (
            *(round(float(v), dedup_decimals) for v in state0),
            round(float(period), dedup_decimals),
            round(float(jacobi), dedup_decimals),
        )
        if row in seen:
            continue
        seen.add(row)
        rows.append(row)
        tags.append(
            heuristic_family_tag(
                state0,
                period,
                jacobi,
                jacobi_bounds=jacobi_bounds,
                period_bounds=period_bounds,
            )
        )
    if not rows:
        raise ValueError(
            "assemble_corpus: no record survived the filters -- check the input "
            "paths and bounds before concluding the corpus is empty."
        )
    features = np.array(rows, dtype=np.float64)
    from cyclerfinder.core.cr3bp import cr3bp_system

    system = cr3bp_system(primary, secondary)
    return OrbitCorpus(
        features=features,
        mu=system.mu,
        primary=primary,
        secondary=secondary,
        n_scanned=n_scanned,
        n_converged_prefilter=n_converged,
        family_tags=np.array(tags, dtype="<U16"),
    )


def _sample_latent_from_clusters(
    cluster_weights: NDArray[np.float64],
    cluster_means: list[NDArray[np.float64]],
    cluster_covs: list[NDArray[np.float64]],
    n_latent: int,
    n: int,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Shared cluster-mixture sampler: draw ``n`` latent points from a
    categorical-over-clusters + per-cluster multivariate-Gaussian mixture.

    Factored out of :meth:`ClusteredGaussianLatentModel.sample` so the SAME
    sampling logic (not a re-derived one) also drives #614's family-
    conditioned and autoencoder-latent variants -- the only thing that
    differs between those variants is how the clusters/latent space were
    built, not how sampling happens once they exist.
    """
    k = len(cluster_weights)
    choices = rng.choice(k, size=n, p=cluster_weights)
    latent = np.empty((n, n_latent))
    for c in range(k):
        idx = np.flatnonzero(choices == c)
        if idx.size == 0:
            continue
        latent[idx] = rng.multivariate_normal(cluster_means[c], cluster_covs[c], size=idx.size)
    return latent


def _fit_cluster_gaussians_from_labels(
    latent: NDArray[np.float64],
    labels: NDArray[np.integer],
    n_clusters: int,
    cov_reg: float,
    centroids: NDArray[np.float64] | None = None,
) -> tuple[list[NDArray[np.float64]], list[NDArray[np.float64]], NDArray[np.float64]]:
    """Shared per-cluster empirical-Gaussian fit given ANY integer partition
    of ``latent`` into ``[0, n_clusters)`` labels (k-means labels for
    :func:`fit_clustered_gaussian`, or heuristic family-tag labels for
    :func:`fit_family_conditioned_gaussian` -- the fitting logic is identical
    either way, only the partition's origin differs).
    """
    n_latent = latent.shape[1]
    identity = np.eye(n_latent) * cov_reg
    global_cov = np.cov(latent.T) + identity
    global_mean = latent.mean(axis=0)
    cluster_means: list[NDArray[np.float64]] = []
    cluster_covs: list[NDArray[np.float64]] = []
    weights: list[float] = []
    for c in range(n_clusters):
        pts = latent[labels == c]
        if len(pts) < n_latent + 2:
            # Too few points in this cluster for a stable covariance estimate;
            # fall back to its centroid (or the global mean if empty) with the
            # GLOBAL covariance -- keeps sampling well-defined without ever
            # collapsing to a zero-variance degenerate cluster.
            if len(pts) == 0:
                fallback_mean = centroids[c] if centroids is not None else global_mean
            else:
                fallback_mean = pts.mean(axis=0)
            cluster_means.append(fallback_mean)
            cluster_covs.append(global_cov)
            weights.append(max(len(pts), 1))
            continue
        cluster_means.append(pts.mean(axis=0))
        cluster_covs.append(np.cov(pts.T) + identity)
        weights.append(len(pts))
    weights_arr = np.asarray(weights, dtype=np.float64)
    weights_arr /= weights_arr.sum()
    return cluster_means, cluster_covs, weights_arr


@dataclass
class ClusteredGaussianLatentModel:
    """Standardize -> linear PCA -> k-means partition -> per-cluster Gaussian.

    The bounded, numpy/scipy-only analog of a CNN-VAE (see module docstring):
    a *linear* encoder/decoder (PCA) instead of a learned nonlinear one, and a
    k-means-partitioned per-cluster empirical Gaussian instead of a single
    learned latent prior. ``sample`` draws NEW feature rows (original units,
    not required to match any single training point) from this fitted
    density; ``transform``/``inverse_transform`` are the encoder/decoder.
    """

    mean: NDArray[np.float64]
    scale: NDArray[np.float64]
    components: NDArray[np.float64]  # (n_features, n_latent)
    explained_variance_ratio: NDArray[np.float64]
    cluster_means: list[NDArray[np.float64]]
    cluster_covs: list[NDArray[np.float64]]
    cluster_weights: NDArray[np.float64]
    n_latent: int
    seed: int = 0
    _rng: np.random.Generator = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._rng = np.random.default_rng(self.seed)

    def transform(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        """Encode original-unit feature rows into latent (PCA) coordinates."""
        return ((features - self.mean) / self.scale) @ self.components

    def inverse_transform(self, latent: NDArray[np.float64]) -> NDArray[np.float64]:
        """Decode latent coordinates back to original-unit feature rows."""
        return (latent @ self.components.T) * self.scale + self.mean

    def sample(self, n: int, *, rng: np.random.Generator | None = None) -> NDArray[np.float64]:
        """Sample ``n`` new feature rows (original units) from the fitted density."""
        rng = rng if rng is not None else self._rng
        latent = _sample_latent_from_clusters(
            self.cluster_weights, self.cluster_means, self.cluster_covs, self.n_latent, n, rng
        )
        return self.inverse_transform(latent)


def fit_clustered_gaussian(
    features: NDArray[np.float64],
    *,
    n_latent: int = 5,
    n_clusters: int = 8,
    cov_reg: float = 1e-3,
    seed: int = 0,
) -> ClusteredGaussianLatentModel:
    """Fit :class:`ClusteredGaussianLatentModel` on ``features`` (n, n_features)."""
    if n_latent > features.shape[1]:
        raise ValueError(f"n_latent={n_latent} exceeds feature dimension {features.shape[1]}")
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale = np.where(scale < 1e-12, 1.0, scale)
    standardized = (features - mean) / scale

    # PCA via SVD -- numpy-only, no sklearn dependency.
    _u, s, vt = np.linalg.svd(standardized, full_matrices=False)
    components = vt[:n_latent].T  # (n_features, n_latent)
    variance = (s**2) / max(len(features) - 1, 1)
    total_variance = variance.sum()
    explained = variance[:n_latent] / total_variance if total_variance > 0 else variance[:n_latent]
    latent = standardized @ components

    centroids, labels = kmeans2(latent, n_clusters, seed=seed, minit="++")
    cluster_means, cluster_covs, weights_arr = _fit_cluster_gaussians_from_labels(
        latent, labels, n_clusters, cov_reg, centroids=centroids
    )

    return ClusteredGaussianLatentModel(
        mean=mean,
        scale=scale,
        components=components,
        explained_variance_ratio=explained,
        cluster_means=cluster_means,
        cluster_covs=cluster_covs,
        cluster_weights=weights_arr,
        n_latent=n_latent,
        seed=seed,
    )


def fit_family_conditioned_gaussian(
    features: NDArray[np.float64],
    family_tags: NDArray[np.str_],
    *,
    n_latent: int = 5,
    cov_reg: float = 1e-3,
    seed: int = 0,
) -> ClusteredGaussianLatentModel:
    """#614: same linear-PCA encoder as :func:`fit_clustered_gaussian`, but

    partition the latent space by the SAME heuristic family tag used to
    condition/stratify (see :func:`heuristic_family_tag`) instead of an
    unsupervised k-means partition. Isolates the "does family-conditioning
    help" question from the "does a nonlinear encoder help" question (this
    still uses the linear PCA encoder) -- the two are independent axes per
    `#614`'s scope. Returns the SAME :class:`ClusteredGaussianLatentModel`
    type as the k-means variant, so both feed the identical downstream
    ``sample``/``transform``/``inverse_transform`` evaluation code.
    """
    if n_latent > features.shape[1]:
        raise ValueError(f"n_latent={n_latent} exceeds feature dimension {features.shape[1]}")
    if len(family_tags) != len(features):
        raise ValueError(
            f"family_tags length {len(family_tags)} does not match features length {len(features)}"
        )
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale = np.where(scale < 1e-12, 1.0, scale)
    standardized = (features - mean) / scale

    _u, s, vt = np.linalg.svd(standardized, full_matrices=False)
    components = vt[:n_latent].T
    variance = (s**2) / max(len(features) - 1, 1)
    total_variance = variance.sum()
    explained = variance[:n_latent] / total_variance if total_variance > 0 else variance[:n_latent]
    latent = standardized @ components

    unique_tags, labels = np.unique(family_tags, return_inverse=True)
    n_clusters = len(unique_tags)
    cluster_means, cluster_covs, weights_arr = _fit_cluster_gaussians_from_labels(
        latent, labels, n_clusters, cov_reg
    )

    return ClusteredGaussianLatentModel(
        mean=mean,
        scale=scale,
        components=components,
        explained_variance_ratio=explained,
        cluster_means=cluster_means,
        cluster_covs=cluster_covs,
        cluster_weights=weights_arr,
        n_latent=n_latent,
        seed=seed,
    )


@dataclass
class AutoencoderModel:
    """Standardize -> shallow tanh-MLP encoder -> linear bottleneck -> shallow

    tanh-MLP decoder. #614's bounded NONLINEAR alternative to
    :class:`ClusteredGaussianLatentModel`'s linear-PCA encoder (`#608`'s own
    flagged limitation #2: a linear encoder may not tightly track the true
    curved family manifolds). Deliberately a plain from-scratch numpy/scipy
    autoencoder fit with ``scipy.optimize.minimize`` (L-BFGS-B, analytic
    gradients) rather than a new torch/jax dependency -- kernel PCA was the
    other option `#614` sanctioned, but this project has no sklearn
    dependency either (checked ``pyproject.toml``) AND kernel PCA has no
    natural decoder (the "pre-image problem"), which this comparison needs
    for generation; a small autoencoder gives both an encoder and an exact
    decoder for free. Same ``transform``/``inverse_transform`` contract as
    the PCA model so it's a drop-in encoder for the SAME clustering +
    sampling machinery.
    """

    mean: NDArray[np.float64]
    scale: NDArray[np.float64]
    w1: NDArray[np.float64]
    b1: NDArray[np.float64]
    w2: NDArray[np.float64]
    b2: NDArray[np.float64]
    w3: NDArray[np.float64]
    b3: NDArray[np.float64]
    w4: NDArray[np.float64]
    b4: NDArray[np.float64]
    n_latent: int
    train_reconstruction_mse: float

    def transform(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        """Encode original-unit feature rows into latent (bottleneck) coordinates."""
        x = (features - self.mean) / self.scale
        h1 = np.tanh(x @ self.w1 + self.b1)
        return h1 @ self.w2 + self.b2

    def inverse_transform(self, latent: NDArray[np.float64]) -> NDArray[np.float64]:
        """Decode latent coordinates back to original-unit feature rows."""
        h2 = np.tanh(latent @ self.w3 + self.b3)
        xhat = h2 @ self.w4 + self.b4
        return xhat * self.scale + self.mean


def fit_autoencoder(
    features: NDArray[np.float64],
    *,
    n_latent: int = 5,
    n_hidden: int = 16,
    l2: float = 1e-4,
    maxiter: int = 300,
    seed: int = 0,
) -> AutoencoderModel:
    """Fit a shallow nonlinear :class:`AutoencoderModel` on ``features`` (n, n_features).

    Architecture: standardize -> tanh(x @ w1 + b1) -> linear bottleneck
    (h1 @ w2 + b2) -> tanh(z @ w3 + b3) -> linear output (h2 @ w4 + b4),
    trained end-to-end by minimizing mean-squared reconstruction error plus
    a small L2 weight penalty via full-batch L-BFGS-B with analytically
    derived (hand-written backprop) gradients -- no autodiff dependency
    needed for a network this small (~n_features*n_hidden*2 + n_hidden*n_latent*2
    weights).
    """
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale = np.where(scale < 1e-12, 1.0, scale)
    x = (features - mean) / scale
    n, n_features = x.shape

    rng = np.random.default_rng(seed)
    shapes: list[tuple[str, tuple[int, ...]]] = [
        ("w1", (n_features, n_hidden)),
        ("b1", (n_hidden,)),
        ("w2", (n_hidden, n_latent)),
        ("b2", (n_latent,)),
        ("w3", (n_latent, n_hidden)),
        ("b3", (n_hidden,)),
        ("w4", (n_hidden, n_features)),
        ("b4", (n_features,)),
    ]
    sizes = [int(np.prod(shape)) for _, shape in shapes]
    split_idx = np.cumsum(sizes)[:-1]

    def unpack(theta: NDArray[np.float64]) -> dict[str, NDArray[np.float64]]:
        parts = np.split(theta, split_idx)
        return {
            name: part.reshape(shape) for (name, shape), part in zip(shapes, parts, strict=True)
        }

    def forward(
        p: dict[str, NDArray[np.float64]],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        h1 = np.tanh(x @ p["w1"] + p["b1"])
        z = h1 @ p["w2"] + p["b2"]
        h2 = np.tanh(z @ p["w3"] + p["b3"])
        xhat = h2 @ p["w4"] + p["b4"]
        return h1, z, h2, xhat

    weight_names = {"w1", "w2", "w3", "w4"}

    def loss_and_grad(theta: NDArray[np.float64]) -> tuple[float, NDArray[np.float64]]:
        p = unpack(theta)
        h1, z, h2, xhat = forward(p)
        diff = xhat - x
        mse = float(np.mean(np.sum(diff**2, axis=1)))
        l2_term = float(l2 * sum(np.sum(p[name] ** 2) for name in weight_names))
        loss = mse + l2_term

        dxhat = (2.0 / n) * diff
        dw4 = h2.T @ dxhat + 2 * l2 * p["w4"]
        db4 = dxhat.sum(axis=0)
        dh2 = dxhat @ p["w4"].T
        dpre2 = dh2 * (1.0 - h2**2)
        dw3 = z.T @ dpre2 + 2 * l2 * p["w3"]
        db3 = dpre2.sum(axis=0)
        dz = dpre2 @ p["w3"].T
        dw2 = h1.T @ dz + 2 * l2 * p["w2"]
        db2 = dz.sum(axis=0)
        dh1 = dz @ p["w2"].T
        dpre1 = dh1 * (1.0 - h1**2)
        dw1 = x.T @ dpre1 + 2 * l2 * p["w1"]
        db1 = dpre1.sum(axis=0)

        grad = np.concatenate([g.ravel() for g in (dw1, db1, dw2, db2, dw3, db3, dw4, db4)])
        return loss, grad

    theta0 = np.concatenate([rng.normal(scale=0.3, size=shape).ravel() for _, shape in shapes])
    result = minimize(
        loss_and_grad, theta0, jac=True, method="L-BFGS-B", options={"maxiter": maxiter}
    )
    params = unpack(result.x)
    _h1, _z, _h2, xhat_final = forward(params)
    final_mse = float(np.mean(np.sum((xhat_final - x) ** 2, axis=1)))

    return AutoencoderModel(
        mean=mean,
        scale=scale,
        w1=params["w1"],
        b1=params["b1"],
        w2=params["w2"],
        b2=params["b2"],
        w3=params["w3"],
        b3=params["b3"],
        w4=params["w4"],
        b4=params["b4"],
        n_latent=n_latent,
        train_reconstruction_mse=final_mse,
    )


@dataclass
class AutoencoderClusteredGaussianModel:
    """#614: k-means-partitioned per-cluster Gaussian over an

    :class:`AutoencoderModel`'s (nonlinear) latent space -- the direct
    nonlinear-encoder analog of :class:`ClusteredGaussianLatentModel`, using
    the SAME cluster-sampling machinery (:func:`_sample_latent_from_clusters`)
    so the only thing that differs from the linear baseline is the
    encoder/decoder.
    """

    encoder: AutoencoderModel
    cluster_means: list[NDArray[np.float64]]
    cluster_covs: list[NDArray[np.float64]]
    cluster_weights: NDArray[np.float64]
    n_latent: int
    seed: int = 0
    _rng: np.random.Generator = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._rng = np.random.default_rng(self.seed)

    def transform(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.encoder.transform(features)

    def inverse_transform(self, latent: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.encoder.inverse_transform(latent)

    def sample(self, n: int, *, rng: np.random.Generator | None = None) -> NDArray[np.float64]:
        rng = rng if rng is not None else self._rng
        latent = _sample_latent_from_clusters(
            self.cluster_weights, self.cluster_means, self.cluster_covs, self.n_latent, n, rng
        )
        return self.inverse_transform(latent)


def fit_autoencoder_clustered_gaussian(
    features: NDArray[np.float64],
    *,
    n_latent: int = 5,
    n_clusters: int = 8,
    n_hidden: int = 16,
    cov_reg: float = 1e-3,
    l2: float = 1e-4,
    maxiter: int = 300,
    seed: int = 0,
) -> AutoencoderClusteredGaussianModel:
    """Fit :class:`AutoencoderClusteredGaussianModel`: nonlinear autoencoder

    encoder + k-means partition + per-cluster Gaussian in its latent space.
    Mirrors :func:`fit_clustered_gaussian`'s structure exactly, substituting
    only the encoder (autoencoder instead of linear PCA) -- so `#614`'s
    linear-vs-nonlinear comparison isolates that one variable.
    """
    encoder = fit_autoencoder(
        features, n_latent=n_latent, n_hidden=n_hidden, l2=l2, maxiter=maxiter, seed=seed
    )
    latent = encoder.transform(features)
    centroids, labels = kmeans2(latent, n_clusters, seed=seed, minit="++")
    cluster_means, cluster_covs, weights_arr = _fit_cluster_gaussians_from_labels(
        latent, labels, n_clusters, cov_reg, centroids=centroids
    )
    return AutoencoderClusteredGaussianModel(
        encoder=encoder,
        cluster_means=cluster_means,
        cluster_covs=cluster_covs,
        cluster_weights=weights_arr,
        n_latent=n_latent,
        seed=seed,
    )


def nearest_neighbor_distances(
    query: NDArray[np.float64], reference: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Euclidean distance from each row of ``query`` to its nearest row in ``reference``.

    Used to check whether a generated/refined candidate is "meaningfully new"
    (far from every training example) or a trivial near-duplicate. Both
    arrays should already be in the SAME normalized (e.g. standardized) space
    for the distances to be comparable across features of different scale.
    """
    tree = cKDTree(reference)
    dist, _ = tree.query(query, k=1)
    return np.asarray(dist, dtype=np.float64)


def uniform_bounding_box_sample(
    features: NDArray[np.float64], n: int, *, rng: np.random.Generator
) -> NDArray[np.float64]:
    """Naive baseline: draw ``n`` points uniformly from ``features``' per-dimension bounding box.

    This is the "no learned structure" comparison point -- same domain the
    generative model was trained on, but with no density information at all.
    """
    lo = features.min(axis=0)
    hi = features.max(axis=0)
    return rng.uniform(lo, hi, size=(n, features.shape[1]))
