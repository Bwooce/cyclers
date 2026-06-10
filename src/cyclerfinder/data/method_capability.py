"""Method-capability partial order + re-sweep gate (Forge Phase 6, design §6a/§6b).

A bare git SHA does NOT express what a method can *reach*: two SHAs of the
single-ellipse genome are equal in capability, while the multi-arc genome at any
SHA reaches strictly more of the trajectory space. So every empty-region record
(``data/empty_regions.jsonl``) carries a :class:`MethodCapability` descriptor —
the genome/corrector plus a capability-tag set — and the re-sweep gate decides,
by capability *subsumption* (not region-match alone), whether a recorded
"empty" still binds a proposed method.

The binding invariant (design §6a): **"region X is empty" always means "empty as
far as method M could reach"** — never an absolute claim.

The partial order (design §6b, ``_CAPABILITY_EDGES``): ``stronger ⊐ weaker``:

* ``multi-arc`` ⊐ ``single-arc`` (more arcs reach more geometry; #163 ⊐ #137),
* ``n-body`` ⊐ ``patched-conic`` (full force supersedes the conic seed),
* ``powered`` / ``low-thrust`` ⊐ ``ballistic`` (added control DOF),
* ``one-dsm-per-leg`` ⊐ ``single-arc`` (added per-leg DOF),
* ``broken-plane`` ⊐ ``coplanar`` (out-of-plane DOF),
* ``leveraging`` ⊐ ``single-arc`` (VILM resonant-leg DOF ⊐ no-leveraging).

NON-GOLDEN / sourced-discipline: the edge set is a *design decision* transcribed
from §6b, not a computed value. It lives in one named constant so the partial
order is auditable in one place.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cyclerfinder.data.empty_regions import EmptyRegionReport


@dataclass(frozen=True)
class MethodCapability:
    """What a search method can *reach* (design §6a) — the empty-region descriptor.

    Attributes
    ----------
    genome:
        The trajectory genome, e.g. ``"single-ellipse free-return"``.
    corrector:
        The corrector/solver, e.g. ``"ballistic_correct"``.
    capability_tags:
        The capability envelope as a tag set (partially ordered by
        :data:`_CAPABILITY_EDGES`).
    git_sha:
        The producing commit (reproducibility; NOT the capability — that is the
        tag set).
    """

    genome: str
    corrector: str
    capability_tags: frozenset[str]
    git_sha: str


# The capability partial order as a "stronger ⊐ weaker" edge set (design §6b).
# Each (a, b) means a is STRICTLY MORE CAPABLE than b: a method carrying tag `a`
# reaches everything tag `b` reaches. Audit the whole order here, in one place.
_CAPABILITY_EDGES: frozenset[tuple[str, str]] = frozenset(
    {
        ("multi-arc", "single-arc"),
        ("n-body", "patched-conic"),
        ("powered", "ballistic"),
        ("low-thrust", "ballistic"),
        ("one-dsm-per-leg", "single-arc"),
        ("broken-plane", "coplanar"),
        ("leveraging", "single-arc"),  # VILM resonant-leg DOF ⊐ no-leveraging
    }
)


def _reaches(tag: str) -> frozenset[str]:
    """The set of tags a method carrying ``tag`` can reach (transitive closure).

    ``tag`` itself plus every tag weaker than it under :data:`_CAPABILITY_EDGES`.
    """
    reached = {tag}
    frontier = [tag]
    while frontier:
        current = frontier.pop()
        for stronger, weaker in _CAPABILITY_EDGES:
            if stronger == current and weaker not in reached:
                reached.add(weaker)
                frontier.append(weaker)
    return frozenset(reached)


def _envelope(tags: Iterable[str]) -> frozenset[str]:
    """The full reachable envelope of a tag set (union of each tag's reach)."""
    env: set[str] = set()
    for tag in tags:
        env |= _reaches(tag)
    return frozenset(env)


def subsumes(a: MethodCapability, b: MethodCapability) -> bool:
    """Does method ``a`` subsume method ``b`` (``a``'s envelope ⊇ ``b``'s)?

    True iff every capability tag of ``b`` is reached by ``a`` under the partial
    order — i.e. ``b``'s envelope is contained in ``a``'s. Reflexive
    (``subsumes(a, a) is True``) and returns ``False`` for incomparable methods
    (neither envelope contains the other).
    """
    a_env = _envelope(a.capability_tags)
    return _envelope(b.capability_tags) <= a_env


def should_sweep(
    *,
    region_id: str,
    method: MethodCapability,
    registry: Sequence[EmptyRegionReport],
) -> bool:
    """The capability-subsumption re-sweep gate (design §6b).

    SKIP the region (return ``False``) ONLY IF a prior empty-region record over
    the same ``region_id`` carries a method-capability that **subsumes** the
    proposed ``method`` — a weaker-or-equal method re-running a region a stronger
    method already emptied learns nothing new.

    RE-SWEEP (return ``True``) when no prior subsumes the proposed method: a
    new/more-capable OR an incomparable method reaches ground the old one could
    not (the #163-reopens-#137 lesson; incomparable never skips), and an unswept
    region always sweeps.
    """
    priors = [r for r in registry if r.region_id == region_id]
    return all(not subsumes(prior.method_capability, method) for prior in priors)


__all__ = ["MethodCapability", "should_sweep", "subsumes"]
