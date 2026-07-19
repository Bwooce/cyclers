"""GAIO-style set-oriented transfer-operator machinery (task #664, `#661` shortlist
item 2).

A genuinely new search paradigm for this project, distinct from every method
built so far (symmetric-closure enumeration, shooting/collocation correctors,
generative-ML seeding, interval certification): instead of searching for
individual periodic orbits or tori, subdivide a region of phase space into
small boxes, estimate a discrete transfer (Perron-Frobenius/Ulam) operator by
short forward propagations of sample points drawn from each box, and identify
ALMOST-INVARIANT SETS -- regions that trap probability mass for a long time
relative to the rest of the covered domain -- from the leading eigenvectors of
the resulting sparse transition matrix. This module is deliberately
model-agnostic: it knows nothing about CR3BP, Poincare sections, or any
specific dynamical system. The Sun-Jupiter Dellnitz-2005 positive control
lives in ``search/quasi_hilda_positive_control.py`` and is built entirely on
top of the generic primitives here.

Method lineage / an honest note on the eigenvector-vs-graph-partition choice
-----------------------------------------------------------------------------
Dellnitz, Junge, Lo, Marsden, Padberg, Preis, Ross & Thiere, "Transport of
Mars-Crossing Asteroids from the Quasi-Hilda Region," Phys. Rev. Lett. 94
(2005) 231102 -- the paper this task's positive control reproduces -- uses
GAIO's box-covering + short-time-propagation transition-matrix machinery
(implemented here, faithfully) but then partitions the resulting transition
GRAPH with the PARTY multilevel graph-partitioning library (a separate,
academic, non-Python dependency) to find almost-invariant sets. This task's
own dispatch instructions direct extraction "from leading eigenvectors" of
the sparse transition matrix instead -- the ORIGINAL GAIO almost-invariant-set
method (Dellnitz & Junge, "On the Approximation of Complicated Dynamical
Behavior," SIAM J. Numer. Anal. 36(2) (1999) 491-515), which the 2005 paper's
own graph-partitioning approach was built to REFINE, not replace. Both
methods are standard GAIO-lineage tools for the same underlying object
(almost-invariant sets of a transfer operator); this module implements the
eigenvector-sign approach as directed, and this docstring flags the
substitution explicitly rather than silently presenting it as identical to
Dellnitz 2005's own exact algorithm.

Escape / leakage convention
----------------------------
The transition matrix built here is generally SUB-stochastic: a sample point
that maps outside the covered box collection (escapes the domain entirely)
contributes no destination-box mass. This is the "open" or "leaky" transfer-
operator convention (Froyland's extension of the closed-system Dellnitz-Junge
formalism) rather than an artificial reflecting boundary -- appropriate here
because box-grid domains are deliberately bounded to a finite phase-space
region and genuine leakage across that boundary is a real, physically
meaningful feature of the dynamics (not a numerical artifact to be papered
over with a self-loop). The leading eigenvalue of a leaky P is < 1 and serves
directly as a retention/escape-rate diagnostic.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray


@dataclass(frozen=True)
class BoxGrid:
    """A regular rectangular grid of boxes covering a ``d``-dimensional box.

    ``lower``/``upper`` are the domain bounds (each shape ``(d,)``);
    ``counts`` is the number of boxes per dimension. Boxes are indexed by a
    flat integer in ``[0, n_boxes)`` via row-major (C order) unraveling of the
    per-dimension multi-index, matching ``numpy.ravel_multi_index`` /
    ``numpy.unravel_index`` conventions directly (used internally).
    """

    lower: NDArray[np.float64]
    upper: NDArray[np.float64]
    counts: tuple[int, ...]

    def __post_init__(self) -> None:
        d = len(self.counts)
        if self.lower.shape != (d,) or self.upper.shape != (d,):
            raise ValueError(
                f"BoxGrid: lower/upper must have shape ({d},) matching len(counts)={d}, "
                f"got lower.shape={self.lower.shape}, upper.shape={self.upper.shape}"
            )
        if np.any(self.upper <= self.lower):
            raise ValueError(
                f"BoxGrid: upper must exceed lower in every dimension, got "
                f"lower={self.lower}, upper={self.upper}"
            )
        if any(c < 1 for c in self.counts):
            raise ValueError(f"BoxGrid: all counts must be >= 1, got {self.counts}")

    @property
    def dim(self) -> int:
        return len(self.counts)

    @property
    def n_boxes(self) -> int:
        n = 1
        for c in self.counts:
            n *= c
        return n

    def widths(self) -> NDArray[np.float64]:
        return (self.upper - self.lower) / np.asarray(self.counts, dtype=np.float64)

    def flat_to_multi(self, flat: int) -> tuple[int, ...]:
        return tuple(int(v) for v in np.unravel_index(flat, self.counts))

    def multi_to_flat(self, multi: tuple[int, ...]) -> int:
        return int(np.ravel_multi_index(multi, self.counts))

    def box_bounds(self, flat: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return ``(lo, hi)`` bounds (each shape ``(d,)``) of box ``flat``."""
        multi = self.flat_to_multi(flat)
        w = self.widths()
        lo = self.lower + np.asarray(multi, dtype=np.float64) * w
        hi = lo + w
        return lo, hi

    def box_center(self, flat: int) -> NDArray[np.float64]:
        lo, hi = self.box_bounds(flat)
        return (lo + hi) / 2.0

    def index_of_point(self, point: NDArray[np.float64]) -> int | None:
        """Flat box index containing ``point``, or ``None`` if outside the grid."""
        p = np.asarray(point, dtype=np.float64)
        if p.shape != (self.dim,):
            raise ValueError(f"BoxGrid.index_of_point: expected shape ({self.dim},), got {p.shape}")
        if np.any(p < self.lower) or np.any(p >= self.upper):
            return None
        w = self.widths()
        multi = np.minimum(np.floor((p - self.lower) / w).astype(int), np.asarray(self.counts) - 1)
        return self.multi_to_flat(tuple(int(v) for v in multi))

    def sample_uniform(self, flat: int, n: int, rng: np.random.Generator) -> NDArray[np.float64]:
        """``n`` points drawn uniformly at random from box ``flat`` (shape ``(n, d)``)."""
        lo, hi = self.box_bounds(flat)
        return rng.uniform(lo, hi, size=(n, self.dim))


@dataclass(frozen=True)
class TransitionMatrixResult:
    """Output of :func:`build_transition_matrix`.

    ``p`` is a (generally SUB-stochastic -- see module docstring) ``b x b``
    sparse row matrix: ``p[i, j]`` is the estimated probability that a point
    drawn uniformly from box ``i`` lands in box ``j`` after one application of
    the map. ``escape_fraction[i]`` is the fraction of box ``i``'s samples
    that left the covered (masked-in) domain entirely; ``n_samples[i]`` is the
    number of samples actually drawn for box ``i`` (0 for masked-out boxes).
    """

    grid: BoxGrid
    mask: NDArray[np.bool_]
    p: sp.csr_matrix
    escape_fraction: NDArray[np.float64]
    n_samples: NDArray[np.intp]
    n_samples_per_box: int


def build_transition_matrix(
    grid: BoxGrid,
    mask: NDArray[np.bool_],
    map_fn: Callable[[NDArray[np.float64]], NDArray[np.float64] | None],
    n_samples_per_box: int,
    rng: np.random.Generator,
    *,
    max_point_retries: int = 10,
) -> TransitionMatrixResult:
    """Build the sparse Ulam/Perron-Frobenius transition matrix.

    For every masked-in box, draw ``n_samples_per_box`` points uniformly at
    random from that box (retrying up to ``max_point_retries`` times per draw
    if ``map_fn`` signals the drawn point itself is invalid by returning
    ``None`` immediately for it -- see below), apply ``map_fn`` once, and bin
    the image into a destination box. A destination that falls outside the
    grid bounds, in a masked-out box, or for which ``map_fn`` itself returns
    ``None`` (the caller's signal that propagation failed or the point
    escaped some larger physical domain the grid doesn't fully capture) all
    count identically as "escaped" mass for that source box -- tracked in
    ``escape_fraction``, never silently renormalized away.

    ``map_fn`` receives a single point (shape ``(d,)``) and must return either
    the image point (shape ``(d,)``) or ``None``. The SAME ``map_fn`` also
    implicitly encodes "this initial point itself was invalid" by returning
    ``None`` -- callers whose domain has invalid interior points (e.g. a
    physical energy manifold that doesn't fill the bounding box exactly)
    should have ``map_fn`` return ``None`` for those too; this function then
    retries a fresh uniform draw from the same box up to ``max_point_retries``
    times before giving up and counting that sample slot as escaped.
    """
    b = grid.n_boxes
    if mask.shape != (b,):
        raise ValueError(f"build_transition_matrix: mask shape {mask.shape} != ({b},)")
    if n_samples_per_box < 1:
        raise ValueError(
            f"build_transition_matrix: n_samples_per_box must be >= 1, got {n_samples_per_box}"
        )

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    escape_fraction = np.zeros(b, dtype=np.float64)
    n_samples = np.zeros(b, dtype=np.intp)

    for i in range(b):
        if not mask[i]:
            continue
        n_samples[i] = n_samples_per_box
        dest_counts: dict[int, int] = {}
        n_escaped = 0
        for _ in range(n_samples_per_box):
            image = None
            for _retry in range(max_point_retries):
                pt = grid.sample_uniform(i, 1, rng)[0]
                image = map_fn(pt)
                if image is not None:
                    break
            if image is None:
                n_escaped += 1
                continue
            j = grid.index_of_point(image)
            if j is None or not mask[j]:
                n_escaped += 1
                continue
            dest_counts[j] = dest_counts.get(j, 0) + 1
        escape_fraction[i] = n_escaped / n_samples_per_box
        for j, count in dest_counts.items():
            rows.append(i)
            cols.append(j)
            data.append(count / n_samples_per_box)

    p = sp.csr_matrix(
        (
            np.asarray(data, dtype=np.float64),
            (np.asarray(rows, dtype=np.intp), np.asarray(cols, dtype=np.intp)),
        ),
        shape=(b, b),
    )
    return TransitionMatrixResult(
        grid=grid,
        mask=mask,
        p=p,
        escape_fraction=escape_fraction,
        n_samples=n_samples,
        n_samples_per_box=n_samples_per_box,
    )


@dataclass(frozen=True)
class AlmostInvariantDecomposition:
    """Output of :func:`almost_invariant_sets_spectral`.

    ``eigenvalues``/``eigenvectors`` are the top ``n_eigvecs`` (by real part)
    LEFT eigenvectors of ``p`` (computed as right eigenvectors of ``p.T`` --
    the correct object for density evolution under the row-stochastic-style
    convention ``rho_{n+1} = rho_n @ P``), sorted by descending real part of
    the eigenvalue; index 0 is the Perron/leading eigenpair (eigenvalue < 1
    when ``p`` is leaky -- see module docstring -- and is the natural escape-
    rate diagnostic; its eigenvector is the quasi-stationary density, and it
    is EXCLUDED from the sign-pattern clustering below since it does not
    discriminate between almost-invariant sets).

    ``labels`` assigns every masked-in box to one of the discovered clusters
    (an integer in ``[0, n_clusters)``) via the sign pattern of eigenvectors
    ``1..n_eigvecs-1`` (Dellnitz & Junge 1999's sign-structure method,
    generalized from a single bipartition to multiple eigenvectors: each
    distinct combination of signs is one candidate cluster). Masked-out boxes
    get label ``-1``. ``cluster_sizes`` gives the box count of the
    ``n_clusters`` largest clusters actually observed, largest first.
    """

    eigenvalues: NDArray[np.complex128]
    eigenvectors: NDArray[np.complex128]
    labels: NDArray[np.intp]
    cluster_sizes: NDArray[np.intp]


def almost_invariant_sets_spectral(
    result: TransitionMatrixResult,
    n_eigvecs: int = 3,
    *,
    dense_threshold: int = 2000,
) -> AlmostInvariantDecomposition:
    """Extract almost-invariant sets from the leading eigenvectors of ``result.p``.

    Uses a dense ``scipy.linalg.eig`` for small problems (``b <=
    dense_threshold``, the common case for a coarse pilot grid, and more
    numerically robust than ARPACK for small/nearly-decomposable matrices)
    and ``scipy.sparse.linalg.eigs`` (largest real part) otherwise. Rows/cols
    of ``p`` restricted to masked-out boxes are all-zero and contribute a
    degenerate zero-eigenvalue subspace; ``n_eigvecs`` is silently capped so
    ARPACK's ``k < N - 1`` requirement is always satisfiable.
    """
    b = result.p.shape[0]
    mask = result.mask
    n_valid = int(mask.sum())
    if n_valid < n_eigvecs + 1:
        raise ValueError(
            f"almost_invariant_sets_spectral: only {n_valid} masked-in boxes, "
            f"need > n_eigvecs={n_eigvecs}"
        )
    pt = result.p.T.tocsc()
    if b <= dense_threshold:
        import scipy.linalg as sla

        eigvals_all, eigvecs_all = sla.eig(pt.toarray())
    else:
        import scipy.sparse.linalg as spla

        k = min(n_eigvecs + 1, b - 2)
        eigvals_all, eigvecs_all = spla.eigs(pt, k=k, which="LR")

    order = np.argsort(-eigvals_all.real)
    eigvals_all = eigvals_all[order]
    eigvecs_all = eigvecs_all[:, order]
    k_take = min(n_eigvecs, eigvals_all.shape[0])
    eigenvalues = eigvals_all[:k_take]
    eigenvectors = eigvecs_all[:, :k_take]

    labels = np.full(b, -1, dtype=np.intp)
    if k_take >= 2:
        sign_vecs = np.sign(eigenvectors[:, 1:k_take].real)
        # Map each distinct sign-pattern row (among masked-in boxes) to a
        # cluster id, largest cluster first.
        valid_idx = np.where(mask)[0]
        patterns = sign_vecs[valid_idx]
        unique_patterns, inverse = np.unique(patterns, axis=0, return_inverse=True)
        counts = np.bincount(inverse.ravel())
        order_by_size = np.argsort(-counts)
        remap = np.empty(len(unique_patterns), dtype=np.intp)
        for new_id, old_id in enumerate(order_by_size):
            remap[old_id] = new_id
        labels[valid_idx] = remap[inverse.ravel()]
        cluster_sizes = counts[order_by_size]
    else:
        labels[mask] = 0
        cluster_sizes = np.array([n_valid], dtype=np.intp)

    return AlmostInvariantDecomposition(
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        labels=labels,
        cluster_sizes=cluster_sizes,
    )


def almost_invariance_ratio(result: TransitionMatrixResult, box_indices: NDArray[np.intp]) -> float:
    """Fraction of a candidate set ``S``'s own mass that stays in ``S`` after
    one iterate, under a UNIFORM box measure (equal-size boxes, as built by
    :class:`BoxGrid`): ``rho(S) = mean_{i in S} sum_{j in S} P[i, j]``. A
    value near 1 indicates a genuinely almost-invariant set; near 0 indicates
    the candidate set is not dynamically coherent (a spurious sign-pattern
    artifact, not a real metastable structure).
    """
    if box_indices.size == 0:
        return 0.0
    sub = result.p[box_indices][:, box_indices]
    row_sums = np.asarray(sub.sum(axis=1)).ravel()
    return float(row_sums.mean())


def transport_probability(
    result: TransitionMatrixResult,
    source_boxes: NDArray[np.intp],
    target_boxes: NDArray[np.intp],
    n_iters: int,
) -> NDArray[np.float64]:
    """Dellnitz et al. 2005 Eq. (transport probability) -- ``p_{R,Q}(n)`` for
    ``n = 0..n_iters``: starting from a UNIFORM probability distribution over
    ``source_boxes`` (mass 1 total), iterate ``rho_{k+1} = rho_k @ P`` (mass
    that escapes the covered domain -- see module docstring -- simply
    vanishes from the accounting, which is correct: escaped mass is
    definitionally not in ``target_boxes``) and record the fraction of mass in
    ``target_boxes`` at every iterate. Returns an array of shape
    ``(n_iters + 1,)``, index 0 being the (generally near-zero, since
    ``source_boxes`` and ``target_boxes`` are constructed disjoint in the
    positive control) initial overlap.
    """
    b = result.p.shape[0]
    rho = np.zeros(b, dtype=np.float64)
    if source_boxes.size == 0:
        raise ValueError("transport_probability: source_boxes is empty")
    rho[source_boxes] = 1.0 / source_boxes.size
    out = np.empty(n_iters + 1, dtype=np.float64)
    out[0] = rho[target_boxes].sum()
    pt = result.p.T.tocsr()  # (rho @ p)_j == (pt @ rho)_j -- avoids a row/2D-array mypy trap
    for k in range(1, n_iters + 1):
        rho = np.asarray(pt @ rho, dtype=np.float64)
        out[k] = rho[target_boxes].sum()
    return out
