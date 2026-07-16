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


@dataclass(frozen=True)
class OrbitCorpus:
    """Assembled, deduplicated ``(state0, period, jacobi)`` records for one CR3BP pair."""

    features: NDArray[np.float64]  # (n, 8): x0,y0,z0,vx0,vy0,vz0,period,jacobi
    mu: float
    primary: str
    secondary: str
    n_scanned: int
    n_converged_prefilter: int

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
    )


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
        k = len(self.cluster_weights)
        choices = rng.choice(k, size=n, p=self.cluster_weights)
        latent = np.empty((n, self.n_latent))
        for c in range(k):
            idx = np.flatnonzero(choices == c)
            if idx.size == 0:
                continue
            latent[idx] = rng.multivariate_normal(
                self.cluster_means[c], self.cluster_covs[c], size=idx.size
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
    identity = np.eye(n_latent) * cov_reg
    global_cov = np.cov(latent.T) + identity
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
            cluster_means.append(centroids[c] if len(pts) == 0 else pts.mean(axis=0))
            cluster_covs.append(global_cov)
            weights.append(max(len(pts), 1))
            continue
        cluster_means.append(pts.mean(axis=0))
        cluster_covs.append(np.cov(pts.T) + identity)
        weights.append(len(pts))
    weights_arr = np.asarray(weights, dtype=np.float64)
    weights_arr /= weights_arr.sum()

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
