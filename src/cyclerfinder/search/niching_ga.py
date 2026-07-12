"""Deterministic-crowding niching GA layer (#581 stage 1).

Implements the Deterministic Crowding niching method exactly as used by
Gurfil & Kasdin (2002), CMAME 191, 5683-5706 (p. 5687 pseudocode; Mahfoud's
method, their ref [22]): individuals are randomly grouped into parent pairs,
each pair generates two children by standard binary-GA crossover + mutation,
and each child competes against the *closer* parent (not the fitter one) —
the winner of each parent-child tournament moves to the next generation. This
is what lets multiple co-existing optima (orbit families) survive one run,
unlike every existing single-optimum ``differential_evolution`` lane in this
package (``search/optimize.py`` etc.). scipy's DE exposes no per-generation
replacement hook, so the mechanism is implemented directly here as a small
generic layer; fitness functions plug in exactly like a scipy DE objective
(vector in, scalar out), and the existing project penalty conventions apply.

Encoding follows the paper's Table 1: binary strings ("String length 32" =
32 bits per free variable), crossover probability 0.999 (single-point, on the
concatenated string), mutation probability 0.001 per bit. Distance ``d(.)``
in the parent-child pairing rule is Euclidean in the bounds-normalized
phenotype space (the paper does not specify; Mahfoud's deterministic crowding
standardly uses phenotypic distance).

Checkpoint/resume is built in so long runs can be executed in bounded
foreground chunks (per project practice).
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

FitnessFn = Callable[[NDArray[np.float64]], float]


@dataclass(frozen=True)
class DeterministicCrowdingConfig:
    """GA constants. Defaults are Gurfil-Kasdin Table 1 (p. 5689)."""

    population_size: int = 200
    generations: int = 400
    bits_per_variable: int = 32
    p_crossover: float = 0.999
    p_mutation: float = 0.001
    seed: int = 0

    def __post_init__(self) -> None:
        if self.population_size % 2:
            raise ValueError("population_size must be even (random pairing)")


@dataclass
class NichingResult:
    """Final population of a deterministic-crowding run (maximization)."""

    phenotypes: NDArray[np.float64]  # (pop, n_var) decoded incl. fixed vars
    fitness: NDArray[np.float64]  # (pop,)
    generations_run: int
    history: list[dict[str, float]] = field(default_factory=list)

    def best(self) -> tuple[NDArray[np.float64], float]:
        i = int(np.argmax(self.fitness))
        return self.phenotypes[i], float(self.fitness[i])


class _Encoding:
    """Binary genome <-> phenotype mapping over a hyper-rectangle.

    Variables whose bounds are degenerate (lo == hi) are held fixed and
    excluded from the genome, mirroring Gurfil-Kasdin's Table 2 sets where
    most state components are pinned to a constant.
    """

    def __init__(self, bounds: Sequence[tuple[float, float]], bits: int) -> None:
        self.bounds = [(float(lo), float(hi)) for lo, hi in bounds]
        for lo, hi in self.bounds:
            if hi < lo:
                raise ValueError(f"invalid bounds ({lo}, {hi})")
        self.bits = bits
        self.free_idx = [i for i, (lo, hi) in enumerate(self.bounds) if hi > lo]
        self.n_var = len(self.bounds)
        self.n_free = len(self.free_idx)
        self.n_bits = self.n_free * bits
        self._weights = 2.0 ** np.arange(bits - 1, -1, -1)
        self._denom = 2.0**bits - 1.0

    def decode_norm(self, genomes: NDArray[np.bool_]) -> NDArray[np.float64]:
        """(pop, n_bits) bits -> (pop, n_free) values in [0, 1]."""
        pop = genomes.shape[0]
        chunks = genomes.reshape(pop, self.n_free, self.bits).astype(float)
        return np.asarray(chunks @ self._weights / self._denom)

    def decode(self, genomes: NDArray[np.bool_]) -> NDArray[np.float64]:
        """(pop, n_bits) bits -> (pop, n_var) full phenotypes."""
        norm = self.decode_norm(genomes)
        out = np.empty((genomes.shape[0], self.n_var))
        for j, (lo, _hi) in enumerate(self.bounds):
            out[:, j] = lo
        for k, j in enumerate(self.free_idx):
            lo, hi = self.bounds[j]
            out[:, j] = lo + norm[:, k] * (hi - lo)
        return out


def _crossover_mutate(
    g1: NDArray[np.bool_],
    g2: NDArray[np.bool_],
    p_cross: float,
    p_mut: float,
    rng: np.random.Generator,
) -> tuple[NDArray[np.bool_], NDArray[np.bool_]]:
    """Single-point crossover (prob ``p_cross``) + per-bit mutation."""
    c1, c2 = g1.copy(), g2.copy()
    n = c1.size
    if n > 1 and rng.random() < p_cross:
        cut = int(rng.integers(1, n))  # cut in [1, n-1]: both sides non-empty
        c1[cut:], c2[cut:] = g2[cut:].copy(), g1[cut:].copy()
    for c in (c1, c2):
        mask = rng.random(n) < p_mut
        c[mask] = ~c[mask]
    return c1, c2


def run_deterministic_crowding(
    fitness_fn: FitnessFn,
    bounds: Sequence[tuple[float, float]],
    config: DeterministicCrowdingConfig | None = None,
    workers: int = 1,
    checkpoint_path: str | Path | None = None,
    max_generations_this_call: int | None = None,
    progress_fn: Callable[[dict[str, float]], None] | None = None,
) -> NichingResult:
    """Run (or resume) a deterministic-crowding GA, maximizing ``fitness_fn``.

    ``checkpoint_path``: npz file written after every generation; if it exists
    it is resumed (population, fitness, RNG state, generation counter).
    ``max_generations_this_call`` bounds the work done in this invocation so a
    long run can be split across bounded foreground calls.
    """
    if config is None:
        config = DeterministicCrowdingConfig()
    enc = _Encoding(bounds, config.bits_per_variable)
    if enc.n_free == 0:
        raise ValueError("all variables are fixed; nothing to optimize")
    rng = np.random.default_rng(config.seed)
    executor = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    try:
        return _run_loop(
            fitness_fn,
            enc,
            config,
            rng,
            executor,
            checkpoint_path,
            max_generations_this_call,
            progress_fn,
        )
    finally:
        if executor is not None:
            executor.shutdown()


def _run_loop(
    fitness_fn: FitnessFn,
    enc: _Encoding,
    config: DeterministicCrowdingConfig,
    rng: np.random.Generator,
    executor: ProcessPoolExecutor | None,
    checkpoint_path: str | Path | None,
    max_generations_this_call: int | None,
    progress_fn: Callable[[dict[str, float]], None] | None,
) -> NichingResult:
    pop = config.population_size

    ckpt = Path(checkpoint_path) if checkpoint_path is not None else None
    history: list[dict[str, float]] = []
    if ckpt is not None and ckpt.exists():
        data = np.load(ckpt, allow_pickle=False)
        genomes = data["genomes"].astype(bool)
        fitness = data["fitness"].astype(float)
        gen_done = int(data["generation"])
        rng.bit_generator.state = json.loads(str(data["rng_state"]))
        history = json.loads(str(data["history"]))
        if genomes.shape != (pop, enc.n_bits):
            raise ValueError("checkpoint shape mismatch with config/bounds")
    else:
        genomes = rng.random((pop, enc.n_bits)) < 0.5
        fitness = _evaluate(fitness_fn, enc, genomes, executor)
        gen_done = 0

    def save(gen: int) -> None:
        if ckpt is None:
            return
        tmp = ckpt.with_suffix(".tmp.npz")
        np.savez(
            tmp,
            genomes=genomes,
            fitness=fitness,
            generation=gen,
            rng_state=json.dumps(rng.bit_generator.state),
            history=json.dumps(history),
        )
        tmp.replace(ckpt)

    if gen_done == 0:
        save(0)

    remaining = config.generations - gen_done
    if max_generations_this_call is not None:
        remaining = min(remaining, max_generations_this_call)

    for _ in range(remaining):
        gen_done += 1
        perm = rng.permutation(pop)
        children = np.empty_like(genomes)
        for k in range(pop // 2):
            i, j = perm[2 * k], perm[2 * k + 1]
            c1, c2 = _crossover_mutate(
                genomes[i], genomes[j], config.p_crossover, config.p_mutation, rng
            )
            children[2 * k] = c1
            children[2 * k + 1] = c2
        child_fitness = _evaluate(fitness_fn, enc, children, executor)

        parent_norm = enc.decode_norm(genomes)
        child_norm = enc.decode_norm(children)
        for k in range(pop // 2):
            i, j = int(perm[2 * k]), int(perm[2 * k + 1])
            a, b = 2 * k, 2 * k + 1
            d_ii = _dist(parent_norm[i], child_norm[a])
            d_jj = _dist(parent_norm[j], child_norm[b])
            d_ij = _dist(parent_norm[i], child_norm[b])
            d_ji = _dist(parent_norm[j], child_norm[a])
            pairs = ((i, a), (j, b)) if d_ii + d_jj <= d_ij + d_ji else ((i, b), (j, a))
            for pi, ci in pairs:
                if child_fitness[ci] >= fitness[pi]:
                    genomes[pi] = children[ci]
                    fitness[pi] = child_fitness[ci]
                    parent_norm[pi] = child_norm[ci]

        stats = {
            "generation": float(gen_done),
            "fitness_mean": float(np.mean(fitness)),
            "fitness_max": float(np.max(fitness)),
            "fitness_std": float(np.std(fitness)),
        }
        history.append(stats)
        save(gen_done)
        if progress_fn is not None:
            progress_fn(stats)

    return NichingResult(
        phenotypes=enc.decode(genomes),
        fitness=fitness.copy(),
        generations_run=gen_done,
        history=history,
    )


def _dist(a: NDArray[np.float64], b: NDArray[np.float64]) -> float:
    return math.sqrt(float(np.sum((a - b) ** 2)))


def _evaluate(
    fitness_fn: FitnessFn,
    enc: _Encoding,
    genomes: NDArray[np.bool_],
    executor: ProcessPoolExecutor | None,
) -> NDArray[np.float64]:
    """Evaluate a genome batch, deduplicating identical genotypes."""
    phen = enc.decode(genomes)
    keys = [np.packbits(g).tobytes() for g in genomes]
    unique: dict[bytes, int] = {}
    for idx, key in enumerate(keys):
        unique.setdefault(key, idx)
    todo = list(unique.values())
    if executor is not None and len(todo) > 1:
        vals = list(executor.map(fitness_fn, [phen[i] for i in todo], chunksize=4))
    else:
        vals = [fitness_fn(phen[i]) for i in todo]
    by_key = {keys[i]: float(v) for i, v in zip(todo, vals, strict=True)}
    return np.array([by_key[k] for k in keys])
