"""Discovery-campaign engine (Track C, #253).

Generalizes the corpus-search daemon (:mod:`scripts.search_campaign_daemon`)
from a CR3BP screening/data-gen loop into a **discovery** pipeline: for a given
pluggable *search target*, enumerate candidate trajectories, close each to a
single canonical periodicity residual, dedup against the sourced catalogue + the
method-versioned negative registry, and **route** every outcome to the unchanged
SILVER->gauntlet artefacts — never to the catalogue.

Pipeline (design ``docs/notes/2026-06-13-discovery-program-spec.md`` Track C):

1. **Enumerate** candidates in a deterministic order (resumability).
2. **Close** each candidate -> canonical residual (km/s).
3. **Dedup** against (a) ``data/catalogue.yaml`` signatures and (b) the
   negative registry ``data/empty_regions.jsonl``; skip the known.
4. **Route**: residual < gate AND novel -> a SILVER record in
   ``data/review_queue.jsonl`` (the gauntlet entry, NOT the catalogue). A swept
   region yielding nothing -> a method-versioned empty-region record in
   ``data/empty_regions.jsonl``.
5. **Resumable / multi-worker**: a gitignored ``out/`` checkpoint records which
   candidate indices are done; ``--worker-id`` shards the deterministic stream so
   N plain processes use N cores without shared-file contention.

DISCIPLINE (golden, non-negotiable): **NO catalogue writeback, ever.** SILVER ->
review queue (a human promotes, never the loop); a clean negative -> empty-region
registry with a :class:`~cyclerfinder.data.method_capability.MethodCapability`
tag ("empty" is always method-conditional). The first wired target is the #254
repeated-moon multi-rev genome over Jupiter; the :class:`SearchTarget` protocol
keeps the multi-arc harness / other genomes hostable later.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

import numpy as np

from cyclerfinder.data.catalog import Catalog, load_catalog
from cyclerfinder.data.empty_regions import (
    EmptyRegionReport,
    append_empty_region,
    load_empty_regions_list,
)
from cyclerfinder.data.method_capability import MethodCapability, should_sweep
from cyclerfinder.data.review_queue import ReviewQueueEntry, append_review_entry

# ---------------------------------------------------------------------------
# Candidate / closure value objects + the pluggable SearchTarget protocol
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Candidate:
    """One enumerated trajectory candidate (target-agnostic envelope).

    ``index`` is the candidate's position in the target's deterministic
    enumeration (the resumability key). ``signature_hash`` is the canonical
    identity used for dedup. ``payload`` carries the target-specific structural
    description (moon sequence, per-leg resonance, ...) opaque to the engine.
    """

    index: int
    signature_hash: str
    sequence: tuple[str, ...]
    primary: str
    payload: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ClosureResult:
    """Outcome of closing one candidate to its canonical periodicity residual.

    ``residual_kms`` is the ONE canonical scalar the engine routes on (the
    worst V_inf-magnitude continuity defect over the cycle, km/s — same quantity
    the #254 corrector minimises). ``converged`` records whether the closure
    routine completed without a numerical failure; a failed close is NOT a
    discovery and is routed as a non-hit.
    """

    converged: bool
    residual_kms: float
    vinf_per_encounter_kms: tuple[float, ...]
    tof_days: tuple[float, ...]
    detail: dict[str, Any] = field(default_factory=dict)


@runtime_checkable
class SearchTarget(Protocol):
    """A pluggable discovery search target (enumerate -> close -> signature).

    The engine knows nothing about a target's genome; it only drives this
    interface, so the multi-arc harness or any future genome can be hosted by
    implementing it. The repeated-moon genome is the first concrete target
    (:class:`RepeatedMoonTarget`).
    """

    @property
    def target_id(self) -> str:
        """Stable identifier (goes into checkpoint + region ids)."""
        ...

    @property
    def primary(self) -> str:
        """Central body of the swept region (catalogue ``primary`` bucket)."""
        ...

    def method_capability(self) -> MethodCapability:
        """The capability descriptor stamped on any empty-region record."""
        ...

    def enumerate_candidates(self) -> Iterator[Candidate]:
        """Yield candidates in a deterministic, resumable order."""
        ...

    def close(self, candidate: Candidate) -> ClosureResult:
        """Close one candidate to its canonical residual (km/s)."""
        ...


# ---------------------------------------------------------------------------
# Catalogue + negative-registry dedup
# ---------------------------------------------------------------------------

VINF_BIN_KMS: float = 0.05
"""V_inf bin width for the moon-cycler dedup signature (matches the catalogue's
spec-16.2 ``VINF_BIN_KMS``; coarse enough that print-precision noise collapses)."""


def _lex_min_rotation(parts: Sequence[str]) -> tuple[str, ...]:
    """Lexicographically-minimal cyclic rotation of a body-name sequence.

    Mirrors :func:`cyclerfinder.data.catalog._lex_min_rotation` but on a name
    tuple (no dash-join round-trip), so a moon sequence and its rotations
    canonicalise identically (Callisto-Ganymede-... is the same cycle regardless
    of which flyby is "first").
    """
    n = len(parts)
    if n <= 1:
        return tuple(parts)
    best = tuple(parts)
    for i in range(1, n):
        rotated = tuple(parts[i:]) + tuple(parts[:i])
        if rotated < best:
            best = rotated
    return best


def _dedupe_closing_body(seq: Sequence[str]) -> tuple[str, ...]:
    """Drop a trailing body equal to the first (open-sequence convention)."""
    t = tuple(seq)
    if len(t) >= 2 and t[0] == t[-1]:
        return t[:-1]
    return t


def moon_cycler_signature_hash(
    *,
    primary: str,
    sequence: Sequence[str],
    vinf_per_encounter_kms: Sequence[float],
) -> str:
    """Canonical identity hash for a repeated-moon cycler (dedup key).

    Reuses the spec-16.2 recipe (rotation-canonical sequence + V_inf multiset
    binned to :data:`VINF_BIN_KMS`, sorted-key ASCII JSON, sha1) but keyed on the
    moon-tour fields the catalogue actually stores for a Jovicentric row
    (``primary`` + ``sequence_canonical`` + ``vinf_kms_at_encounters``). This is
    the SAME identity used to dedup a closed candidate against the catalogue and
    against itself, so one cycle can never hash two ways.
    """
    seq_canon = _lex_min_rotation(_dedupe_closing_body(sequence))
    vinf_binned = sorted(
        round(float(v) / VINF_BIN_KMS) * VINF_BIN_KMS for v in vinf_per_encounter_kms
    )
    hash_input = {
        "primary": primary,
        "sequence_canonical": "-".join(seq_canon),
        "vinf_multiset_binned": vinf_binned,
    }
    blob = json.dumps(hash_input, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return "sha1:" + hashlib.sha1(blob.encode("utf-8")).hexdigest()


def catalogue_moon_signatures(catalog: Catalog, *, primary: str) -> set[str]:
    """Every moon-cycler dedup signature already in the catalogue for ``primary``.

    Walks the loaded catalogue, keeps the rows whose ``primary`` matches (the
    Jovicentric bucket), and recomputes :func:`moon_cycler_signature_hash` from
    each row's published sequence + per-encounter V_inf. A candidate whose hash
    lands in this set is a known (sourced) cycler -> skipped, never re-queued.
    """
    sigs: set[str] = set()
    for entry in catalog.entries:
        if (entry.primary or "Sun") != primary:
            continue
        vinfs = [v for _, v in entry.vinf_kms_at_encounters if v is not None]
        seq = (
            tuple(entry.sequence_canonical.split("-")) if entry.sequence_canonical else entry.bodies
        )
        if not vinfs or not seq:
            continue
        sigs.add(
            moon_cycler_signature_hash(primary=primary, sequence=seq, vinf_per_encounter_kms=vinfs)
        )
    return sigs


# ---------------------------------------------------------------------------
# The repeated-moon multi-rev search target (the #254 genome, generalized)
# ---------------------------------------------------------------------------

DAY_S: float = 86400.0


def _mean_motion_rad_day(system_mu: float, sma_km: float) -> float:
    """Mean motion (rad/day) about the primary, Kepler III (registry-derived)."""
    return math.sqrt(system_mu / sma_km**3) * DAY_S


def _moon_state(
    theta0: float, n_rad_day: float, t_days: float, sma_km: float, mu: float
) -> tuple[np.ndarray, np.ndarray]:
    """Circular-coplanar moon position/velocity at ``t_days`` (planet frame)."""
    theta = theta0 + n_rad_day * t_days
    v_circ = math.sqrt(mu / sma_km)
    pos = np.array([sma_km * math.cos(theta), sma_km * math.sin(theta), 0.0])
    vel = np.array([-v_circ * math.sin(theta), v_circ * math.cos(theta), 0.0])
    return pos, vel


@dataclass(frozen=True)
class RepeatedMoonTarget:
    """The #254 repeated-moon multi-rev genome as a discovery search target.

    Enumerates bounded moon-sequence x per-leg ``(n_rev)`` combos over a planet's
    registered moons and closes each by phasing the circular-coplanar moons and
    solving the planet-frame Lambert legs (the same machinery the #254 corrector
    / cge_scaffold reproduction use, generalized to arbitrary sequences and
    body-agnostic registry-derived mean motions). The canonical residual is the
    worst per-flyby V_inf-magnitude continuity defect over the cycle (km/s).

    The enumeration is deliberately BOUNDED + deterministic (sorted moon set,
    fixed resonance grid, capped sequence length) so the daemon is resumable and
    the smoke is tiny. The open-ended widening of the grids is the coordinator's
    background run, not the build.
    """

    primary: str = "Jupiter"
    moons: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto")
    seq_lengths: tuple[int, ...] = (3,)
    n_rev_grid: tuple[int, ...] = (0, 1, 2)
    n_phase_samples: int = 12
    tof_resonance_grid: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0)
    git_sha: str = "uncommitted"

    @property
    def target_id(self) -> str:
        return f"repeated-moon-{self.primary.lower()}"

    def method_capability(self) -> MethodCapability:
        """Repeated-moon multi-rev: multi-arc, patched-conic, ballistic, coplanar.

        Tags are the capability envelope (NOT a SHA) consumed by the
        re-sweep gate: this genome reaches multi-arc / multi-rev repeated-moon
        topologies the prior single-ellipse zero-rev genome could not.
        """
        return MethodCapability(
            genome="repeated-moon multi-rev (#254)",
            corrector="periodicity-continuity (planet-frame Lambert)",
            capability_tags=frozenset(
                {"multi-arc", "patched-conic", "ballistic", "coplanar", "leveraging"}
            ),
            git_sha=self.git_sha,
        )

    # -- registry-derived per-moon constants (cached at construction time) --
    def _moon_consts(self) -> dict[str, tuple[float, float]]:
        """``{moon: (sma_km, mean_motion_rad_day)}`` from the satellites registry."""
        from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

        mu = PRIMARIES[self.primary]
        out: dict[str, tuple[float, float]] = {}
        for m in self.moons:
            sat = SATELLITES[m]
            out[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))
        return out

    def _sequences(self) -> list[tuple[str, ...]]:
        """Deterministic bounded list of repeated-moon sequences.

        For each requested length k, every length-k product over the SORTED moon
        set whose consecutive bodies differ (a leg must change moons) and which
        uses at least two distinct moons (a single-moon resonance loop is the
        zero-rev genome's territory, not a repeated-moon tour). Deterministic
        order: ``itertools.product`` over the sorted moons.
        """
        moons = tuple(sorted(self.moons))
        seqs: list[tuple[str, ...]] = []
        for k in self.seq_lengths:
            for combo in itertools.product(moons, repeat=k):
                if any(combo[i] == combo[i + 1] for i in range(k - 1)):
                    continue
                if len(set(combo)) < 2:
                    continue
                seqs.append(combo)
        return seqs

    def enumerate_candidates(self) -> Iterator[Candidate]:
        """Yield candidates in deterministic order (sequence x per-leg n_rev).

        The signature hash at enumeration time uses a *structural* placeholder
        V_inf multiset (zeros) so the index<->candidate map is stable before the
        close; the engine recomputes the true dedup signature from the closed
        V_inf after the close (the structural hash is only the enumeration key).
        """
        idx = 0
        for seq in self._sequences():
            n_legs = len(seq) - 1
            for nrevs in itertools.product(self.n_rev_grid, repeat=n_legs):
                struct_blob = json.dumps(
                    {"seq": list(seq), "n_rev": list(nrevs)},
                    sort_keys=True,
                    separators=(",", ":"),
                )
                struct_hash = "sha1:" + hashlib.sha1(struct_blob.encode("utf-8")).hexdigest()
                yield Candidate(
                    index=idx,
                    signature_hash=struct_hash,
                    sequence=seq,
                    primary=self.primary,
                    payload={"n_rev": list(nrevs)},
                )
                idx += 1

    def close(self, candidate: Candidate) -> ClosureResult:
        """Close one candidate to its canonical V_inf-continuity residual.

        Sweep the bounded phase grid (initial moon longitudes) x the per-leg ToF
        resonance grid; for each phasing, place the circular-coplanar moons at
        the cumulative flyby epochs, solve the planet-frame multi-rev Lambert leg
        at the requested ``n_rev``, and measure the worst per-flyby V_inf-
        magnitude continuity. Return the phasing that MINIMISES that worst defect
        (the corrector's canonical residual). A leg with no Lambert solution at
        the requested revolution count marks the phasing infeasible; if no
        phasing is feasible the candidate is reported non-converged.
        """
        from cyclerfinder.core.lambert import lambert
        from cyclerfinder.core.satellites import PRIMARIES

        mu = PRIMARIES[self.primary]
        consts = self._moon_consts()
        seq = candidate.sequence
        nrevs: list[int] = candidate.payload["n_rev"]
        n_legs = len(seq) - 1

        best_residual = math.inf
        best_vinf: tuple[float, ...] = ()
        best_tof: tuple[float, ...] = ()

        phase0_grid = [
            2.0 * math.pi * i / self.n_phase_samples for i in range(self.n_phase_samples)
        ]
        for phase0 in phase0_grid:
            # Each moon gets a distinct deterministic initial longitude offset.
            theta0 = {
                m: phase0 + 2.0 * math.pi * j / max(1, len(consts))
                for j, m in enumerate(sorted(consts))
            }
            for tof_scale in self.tof_resonance_grid:
                ok, residual, vinfs, tofs = self._close_one_phasing(
                    seq, nrevs, n_legs, theta0, tof_scale, consts, mu, lambert
                )
                if ok and residual < best_residual:
                    best_residual = residual
                    best_vinf = vinfs
                    best_tof = tofs

        if not math.isfinite(best_residual):
            return ClosureResult(
                converged=False,
                residual_kms=math.inf,
                vinf_per_encounter_kms=(),
                tof_days=(),
                detail={"reason": "no feasible phasing (no Lambert solution at requested n_rev)"},
            )
        return ClosureResult(
            converged=True,
            residual_kms=best_residual,
            vinf_per_encounter_kms=best_vinf,
            tof_days=best_tof,
        )

    def _close_one_phasing(
        self,
        seq: tuple[str, ...],
        nrevs: list[int],
        n_legs: int,
        theta0: dict[str, float],
        tof_scale: float,
        consts: dict[str, tuple[float, float]],
        mu: float,
        lambert: Any,
    ) -> tuple[bool, float, tuple[float, ...], tuple[float, ...]]:
        """Close a single (phasing, ToF-scale) sample; return continuity defect.

        ToF per leg is ``tof_scale`` x the geometric-mean orbital period of the
        leg's two moons (a resonance-scaled transfer time bracketing the
        commensurable arc). Returns ``(feasible, worst_continuity_kms,
        vinf_per_flyby, tof_per_leg)``.
        """
        # Per-leg ToF on the resonance grid; cumulative epochs for moon phasing.
        tofs: list[float] = []
        for k in range(n_legs):
            _, na = consts[seq[k]]
            _, nb = consts[seq[k + 1]]
            pa = 2.0 * math.pi / na
            pb = 2.0 * math.pi / nb
            tofs.append(tof_scale * math.sqrt(pa * pb))
        epochs = [0.0]
        for tof in tofs:
            epochs.append(epochs[-1] + tof)

        states = []
        for m, t in zip(seq, epochs, strict=True):
            sma, n = consts[m]
            states.append(_moon_state(theta0[m], n, t, sma, mu))

        vinf_in: list[float | None] = [None] * len(seq)
        vinf_out: list[float | None] = [None] * len(seq)
        for k in range(n_legs):
            r_a, v_a = states[k]
            r_b, v_b = states[k + 1]
            sols = lambert(r_a, r_b, tofs[k] * DAY_S, mu=mu, max_revs=max(0, nrevs[k]))
            wanted = [s for s in sols if s.n_revs == nrevs[k]]
            if not wanted:
                return (False, math.inf, (), ())
            # Pick the lowest-energy (smallest departure V_inf) branch.
            best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
            vinf_out[k] = float(np.linalg.norm(best.v1 - v_a))
            vinf_in[k + 1] = float(np.linalg.norm(best.v2 - v_b))

        # Continuity defect at every interior flyby (both an in and an out side).
        worst = 0.0
        per_flyby: list[float] = []
        for k in range(len(seq)):
            vi = vinf_in[k]
            vo = vinf_out[k]
            if vi is not None and vo is not None:
                worst = max(worst, abs(vi - vo))
            # representative per-flyby V_inf (prefer inbound, else outbound)
            rep = vi if vi is not None else vo
            per_flyby.append(rep if rep is not None else 0.0)
        return (True, worst, tuple(per_flyby), tuple(tofs))


# ---------------------------------------------------------------------------
# The campaign engine: enumerate -> dedup -> close -> route (resumable)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CampaignConfig:
    """Static configuration for one campaign run."""

    gate_residual_kms: float = 0.05
    """Closure must beat this canonical residual to count as a hit (km/s). The
    coordinator tightens it; the build default is the catalogue's V_inf print
    bin so a smoke hit means "closed to within a bin", not a fabricated 1e-11."""

    worker_id: int = 0
    n_workers: int = 1
    max_candidates: int | None = None
    """Hard cap on candidates evaluated this run (the smoke sets it tiny; the
    real run leaves it ``None`` = exhaust the deterministic enumeration)."""


@dataclass
class CampaignRouting:
    """Where outcomes are written (defaults to the real registries).

    The smoke overrides these with temp paths so it NEVER touches the real
    ``data/`` artefacts; ``catalogue`` is read-only and never written.
    """

    review_queue_path: Path
    empty_regions_path: Path
    checkpoint_path: Path


@dataclass
class CampaignStats:
    """Real, reportable counters for a run (no fabrication)."""

    enumerated: int = 0
    evaluated: int = 0
    skipped_done: int = 0
    skipped_known: int = 0
    closed: int = 0
    failed_close: int = 0
    silver_routed: int = 0
    empty_routed: int = 0
    hits: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "enumerated": self.enumerated,
            "evaluated": self.evaluated,
            "skipped_done": self.skipped_done,
            "skipped_known": self.skipped_known,
            "closed": self.closed,
            "failed_close": self.failed_close,
            "silver_routed": self.silver_routed,
            "empty_routed": self.empty_routed,
            "hits": list(self.hits),
        }


def _load_checkpoint(path: Path) -> set[int]:
    """Load the set of completed candidate indices (gitignored ``out/``)."""
    if not path.exists():
        return set()
    done: set[int] = set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line:
            done.add(int(line))
    return done


def _append_checkpoint(path: Path, index: int) -> None:
    """Append one completed candidate index (atomic-enough single-line append)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(f"{index}\n")


def _silver_entry(
    target: SearchTarget,
    candidate: Candidate,
    closure: ClosureResult,
    dedup_hash: str,
) -> ReviewQueueEntry:
    """Build the SILVER review-queue entry for a novel closed candidate.

    The candidate is machine-confirmed (closed below the gate) but UNSOURCED, so
    it is SILVER and lands in the human-review queue with no panel refutation and
    a literature_check left ``None`` (a human checks the literature before any
    promotion). NEVER written to the catalogue.
    """
    return ReviewQueueEntry(
        candidate_id=f"{target.target_id}-{candidate.index:08d}",
        signature_hash=dedup_hash,
        verdict_tier="silver",
        match_outcome="novel",
        known_id=None,
        superseded_by=(),
        vinf_per_encounter_kms=tuple(closure.vinf_per_encounter_kms),
        tof_days=tuple(closure.tof_days),
        bend_feasible=True,
        max_vinf_kms=max(closure.vinf_per_encounter_kms) if closure.vinf_per_encounter_kms else 0.0,
        sequence=tuple(candidate.sequence),
        period_k=len(candidate.sequence) - 1,
        model_assumption="circular-coplanar",
        verdict_audit={
            "residual_kms": closure.residual_kms,
            "n_rev": candidate.payload.get("n_rev"),
            "primary": candidate.primary,
            "closed_by": target.target_id,
        },
        panel={"majority_refute": False, "n_refuted": 0, "n_verifiers": 0},
        t_added="",
        literature_check=None,
    )


def _empty_region_report(
    target: SearchTarget,
    config: CampaignConfig,
    stats: CampaignStats,
) -> EmptyRegionReport:
    """Build the method-versioned empty-region record for a no-hit sweep.

    A swept region yielding zero SILVER survivors is a first-class negative: it
    carries the bounded search extent (candidates evaluated), the prune gates
    (the closure gate + the dedup skips), and the
    :class:`MethodCapability` tag set so a later more-capable method knows to
    re-sweep ("empty" is never unconditional).
    """
    region_id = f"{target.target_id}-sweep"
    return EmptyRegionReport(
        region_id=region_id,
        family="repeated-moon multi-rev cyclers",
        centre=target.primary,
        topologies=(),
        method_capability=target.method_capability(),
        search_extent={
            "points_total": max(1, stats.evaluated),
            "candidates_enumerated": stats.enumerated,
            "candidates_evaluated": stats.evaluated,
            "worker_id": config.worker_id,
            "n_workers": config.n_workers,
        },
        prune_gates=(
            f"closure residual >= gate ({config.gate_residual_kms} km/s)",
            "dedup: known catalogue signature",
            "dedup: prior empty-region (capability-subsumed)",
            "closure failed (no Lambert solution at requested n_rev)",
        ),
        result={
            "silver_routed": stats.silver_routed,
            "closed": stats.closed,
            "failed_close": stats.failed_close,
            "skipped_known": stats.skipped_known,
        },
        verdict="empty",
        interpretation=(
            "No novel repeated-moon cycler closed below the gate in the swept "
            "bounded enumeration; empty as far as this genome could reach."
        ),
        source_anchors="discovery campaign #253 (Track C); no published anchor (novel search)",
        run={
            "target": target.target_id,
            "git_sha": target.method_capability().git_sha,
        },
    )


def run_campaign(
    target: SearchTarget,
    config: CampaignConfig,
    routing: CampaignRouting,
    *,
    catalog: Catalog | None = None,
    empty_registry: Iterable[EmptyRegionReport] | None = None,
    write_empty_on_no_hits: bool = True,
) -> CampaignStats:
    """Run the discovery pipeline for one target (synchronous, resumable).

    enumerate -> (skip done / skip known / skip capability-subsumed) -> close ->
    route. SILVER hits append to ``routing.review_queue_path``; a no-hit sweep
    appends one method-versioned record to ``routing.empty_regions_path``. The
    catalogue is read-only. Progress is checkpointed per candidate to the
    gitignored ``routing.checkpoint_path`` so a restart resumes mid-stream.

    Multi-worker: candidate ``index % n_workers == worker_id`` selects this
    worker's deterministic shard. ``config.max_candidates`` caps the run (the
    smoke uses a tiny cap; the real run leaves it ``None``).
    """
    cat = catalog if catalog is not None else load_catalog()
    known_sigs = catalogue_moon_signatures(cat, primary=target.primary)
    registry = (
        list(empty_registry)
        if empty_registry is not None
        else load_empty_regions_list(routing.empty_regions_path)
    )

    # Capability-subsumption gate: if a prior empty-region over this region was
    # produced by a method that subsumes ours, skip the whole sweep (learns
    # nothing new). A new/incomparable method re-sweeps (the #163-reopens-#137
    # lesson). This is the region-level gate; per-candidate dedup is below.
    sweep_region_id = f"{target.target_id}-sweep"
    region_open = should_sweep(
        region_id=sweep_region_id,
        method=target.method_capability(),
        registry=registry,
    )

    stats = CampaignStats()
    if not region_open:
        return stats

    done = _load_checkpoint(routing.checkpoint_path)

    for candidate in target.enumerate_candidates():
        stats.enumerated += 1
        if candidate.index % config.n_workers != config.worker_id:
            continue
        if candidate.index in done:
            stats.skipped_done += 1
            continue
        if config.max_candidates is not None and stats.evaluated >= config.max_candidates:
            break

        stats.evaluated += 1
        closure = target.close(candidate)
        if not closure.converged:
            stats.failed_close += 1
            _append_checkpoint(routing.checkpoint_path, candidate.index)
            continue
        stats.closed += 1

        # True dedup signature from the CLOSED V_inf (the structural enumeration
        # hash was only the resumability key).
        dedup_hash = moon_cycler_signature_hash(
            primary=candidate.primary,
            sequence=candidate.sequence,
            vinf_per_encounter_kms=closure.vinf_per_encounter_kms,
        )
        if dedup_hash in known_sigs:
            stats.skipped_known += 1
            _append_checkpoint(routing.checkpoint_path, candidate.index)
            continue

        if closure.residual_kms < config.gate_residual_kms:
            entry = _silver_entry(target, candidate, closure, dedup_hash)
            append_review_entry(routing.review_queue_path, entry)
            stats.silver_routed += 1
            stats.hits.append(entry.candidate_id)
            known_sigs.add(dedup_hash)  # self-dedup within the run

        _append_checkpoint(routing.checkpoint_path, candidate.index)

    if stats.silver_routed == 0 and write_empty_on_no_hits and stats.evaluated > 0:
        append_empty_region(routing.empty_regions_path, _empty_region_report(target, config, stats))
        stats.empty_routed += 1

    return stats


__all__ = [
    "CampaignConfig",
    "CampaignRouting",
    "CampaignStats",
    "Candidate",
    "ClosureResult",
    "RepeatedMoonTarget",
    "SearchTarget",
    "catalogue_moon_signatures",
    "moon_cycler_signature_hash",
    "run_campaign",
]
