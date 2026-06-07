"""Multiple-shooting differential corrector in restricted n-body (design §3, #133).

Consumer 2 of the harness: the **Jones flyby-propagation shooter**. Multiple
shooting over encounter nodes — full-state continuity in real (rails) dynamics
between nodes propagated by :class:`cyclerfinder.nbody.propagator.RestrictedNBody`
— that corrects a patched-conic seed to an n-body-ballistic cycler. This is the
SNOPT analogue of Jones, Hernandez & Jesick (AAS 17-577, 2017): their stage 1
accepted conic v∞-mismatches <= 200 m/s, then ran SNOPT n-body correction to
ballistic; here the conic seed is the START and the multiple-shooting solve drives
the full-state defects to zero.

**SEEDING — the #135 verdict (binding).** The #135 like-for-like diagnostic
(``docs/notes/2026-06-06-russell12-likeforlike.md``) found 0 CLOSE-AND-MATCH
coplanar-vs-coplanar on known-solvable instances: the corrector closes
geometrically but lands OFF-ANCHOR (our V∞ 9-28 vs sourced 3-10 km/s). The
verdict is **seeding/basin, not solver deficiency**. The implication, carried into
this shooter, is absolute: a single-shoot from a naive equispaced / coplanar /
blind-scan seed falls into the same high-V∞ basin. **This shooter MUST be seeded
from the #133 near-miss conic survey** — the lowest-residual low-V∞ conic chains
near the Jones anchors collected by :func:`near_miss_survey` — NEVER from a blind
scan. The near-miss survey (Phase 3a) is the seeding source; the multiple-shooting
solve refines from there.

**GOLDEN DISCIPLINE.** The shooter OUTPUT (n-body-ballistic V∞ at the converged
nodes) is the side under test; the EXPECTED side is the SOURCED Jones multiset
(AAS 17-577 Tables 2/3) only. No self-computed value is ever EXPECTED. Divergence
is a first-class non-converged record (mirror ``correct.py``'s honest
non-converged result), never an exception.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

Vec3 = NDArray[np.float64]

# Full Cartesian state dimension per node (r[3] + v[3]).
_STATE_DIM = 6


def defect_count(*, n_encounters: int) -> int:
    """Number of scalar full-state continuity defects for ``n_encounters`` nodes.

    A cycler with ``n`` encounter nodes has ``n - 1`` interior legs, each
    contributing a full 6-component state defect (the propagated arc end vs the
    next node state). E.g. E-M-E-V-V-E (6 encounters) -> 5 legs -> 30 defects.
    """
    if n_encounters < 2:
        raise ValueError("a cycler needs at least 2 encounter nodes")
    return (n_encounters - 1) * _STATE_DIM


def build_shooting_vector(
    nodes: Mapping[str, Vec3],
    epochs: Sequence[float],
    tofs: Sequence[float],
    *,
    slack_leg: int,
    period_days: float,
) -> NDArray[np.float64]:
    """Pack the multiple-shooting free-variable vector ``x``.

    Layout: ``[ {node states, 6 each, in key order}, {node epochs}, {free ToFs} ]``.
    The **slack leg** ToF is eliminated from the free vector (it is reconstructed
    from the period pin ``period_days - sum(other ToFs)``, the
    ``correct.py:_reconstruct_tofs`` convention, ``correct.py:77-86``) so the
    period is held exactly. Node states are stacked in ``sorted`` key order for a
    stable, reproducible packing (the same ``b{i}`` vocabulary as
    ``correct._vinf_nodes``).
    """
    if not 0 <= slack_leg < len(tofs):
        raise ValueError(f"slack_leg {slack_leg} out of range for {len(tofs)} ToFs")
    parts: list[NDArray[np.float64]] = []
    for key in sorted(nodes):
        parts.append(np.asarray(nodes[key], dtype=np.float64).ravel())
    parts.append(np.asarray(epochs, dtype=np.float64).ravel())
    free_tofs = [t for i, t in enumerate(tofs) if i != slack_leg]
    parts.append(np.asarray(free_tofs, dtype=np.float64))
    return np.concatenate(parts)


__all__ = ["build_shooting_vector", "defect_count"]
