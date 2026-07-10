"""Guard against mistaking a numerical artifact for real physics in a sweep.

Found 2026-07-11 (task #559, a follow-up diagnostic on the daily
launch-epoch V4-strict sweep on catalogue row #312): a parameter sweep over
a smooth physical system produced 55 isolated single-point FAIL "spikes"
out of 731 points, each surrounded by PASS on both immediate neighbors. A
first pass reported this as "chaotic/stochastic" without verifying the
cause. Real physical sensitivity in a smooth dynamical system should vary
continuously with small parameter perturbations; an ISOLATED singleton flip
(as opposed to a cluster, a smooth transition band, or a genuine sharp-but-
continuous boundary) is a strong prior for a numerical artifact -- a solver
branch switch, an ephemeris/kernel interpolation-node boundary, a discrete
tolerance edge case -- not genuine chaos.

This module is the reusable, automatic version of "notice the singleton
pattern before trusting it": any future sweep-analysis script (in the
spirit of ``scripts/analyze_338_boundary.py``) should call
:func:`detect_isolated_singleton_anomalies` on its boolean result series and
report the count/fraction explicitly, rather than only reporting an
aggregate pass rate that silently absorbs this pattern.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class SingletonAnomaly:
    """One isolated singleton flip: index ``i`` disagrees with both
    ``i - 1`` and ``i + 1``, which agree with each other."""

    index: int
    value: bool
    label: str | None = None


def detect_isolated_singleton_anomalies(
    values: Sequence[bool],
    labels: Sequence[str] | None = None,
) -> list[SingletonAnomaly]:
    """Return every isolated singleton flip in ``values``.

    A singleton at index ``i`` (``1 <= i <= len(values) - 2``) is isolated
    iff ``values[i] != values[i - 1]``, ``values[i] != values[i + 1]``, and
    ``values[i - 1] == values[i + 1]`` -- i.e. one point disagrees with an
    otherwise-agreeing pair of neighbors on both sides. Endpoints (index 0
    and the last index) are never flagged since they have no neighbor on
    one side to compare against.

    ``labels`` (e.g. epoch strings), if given, must be the same length as
    ``values`` and are carried into the returned :class:`SingletonAnomaly`
    rows for reporting; omit to get bare indices.

    This is a NECESSARY-NOT-SUFFICIENT flag, not a verdict: a real narrow
    dynamical feature (e.g. a genuine one-point-wide resonance crossing at
    exactly this sampling resolution) can also produce an isolated flip.
    Treat every returned anomaly as "verify before trusting," per
    ``[[feedback_isolated_sweep_flips_suspect_artifact]]`` -- zoom into
    finer resolution around it, or trace the code path for discrete/
    branching logic near the anomalous point, before accepting either a
    "real chaotic effect" or "artifact" conclusion.

    CAUTION -- do not concatenate non-contiguous sweep windows before
    calling this. ``values``/``labels`` are treated as genuinely regularly-
    and contiguously-sampled; if two disjoint ranges (e.g. daily epochs
    across 2000, then a SEPARATE daily sweep across 2030) are concatenated
    into one sequence, the boundary between them is treated as if it were a
    real adjacent-sample pair, which can silently manufacture a spurious
    anomaly right at the seam (confirmed empirically in the #559 data: the
    last 2000 row and first 2030 row produced exactly one such artifact).
    Call this once per genuinely contiguous window and sum/report the
    results separately.
    """
    if labels is not None and len(labels) != len(values):
        raise ValueError(f"labels length {len(labels)} != values length {len(values)}")
    anomalies: list[SingletonAnomaly] = []
    for i in range(1, len(values) - 1):
        if (
            values[i] != values[i - 1]
            and values[i] != values[i + 1]
            and values[i - 1] == values[i + 1]
        ):
            label = labels[i] if labels is not None else None
            anomalies.append(SingletonAnomaly(index=i, value=values[i], label=label))
    return anomalies


def singleton_anomaly_summary(
    values: Sequence[bool],
    labels: Sequence[str] | None = None,
) -> str:
    """One-line human-readable summary, suitable for a sweep script's own
    stdout report -- makes the count/fraction impossible to silently miss.
    """
    anomalies = detect_isolated_singleton_anomalies(values, labels)
    n = len(values)
    frac = len(anomalies) / n if n else 0.0
    if not anomalies:
        return f"0/{n} isolated singleton flips -- no artifact-suspicion flag raised."
    sample = ", ".join((a.label if a.label is not None else str(a.index)) for a in anomalies[:5])
    more = f" (+{len(anomalies) - 5} more)" if len(anomalies) > 5 else ""
    return (
        f"{len(anomalies)}/{n} ({frac:.1%}) isolated singleton flips detected -- "
        f"SUSPECT NUMERICAL ARTIFACT, verify before trusting as real physics "
        f"(see feedback_isolated_sweep_flips_suspect_artifact). Examples: {sample}{more}."
    )


__all__ = [
    "SingletonAnomaly",
    "detect_isolated_singleton_anomalies",
    "singleton_anomaly_summary",
]
