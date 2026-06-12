"""Positive-control convention for negative campaigns (FALSE-CONSENSUS DEFENCE, #202).

Doctrine item 4 (``docs/notes/2026-06-11-project-review-results.md`` §"The
false-consensus doctrine"): **each gate must be demonstrated to catch a planted
defect of the class it claims to catch** — and, dually, before a *negative*
result ("method X found nothing / refutes the family") is trusted, the method
must first re-find a KNOWN solution *through the identical pipeline
configuration*. This is the rule that would have killed #180: a method that
cannot rediscover a control it should find is mis-configured, so its empty
result is uninformative, not a real negative.

This module provides the reusable assert helper so a negative-campaign does not
have to hand-roll (and silently weaken) that check. It does NOT retrofit every
campaign; it is opt-in plumbing a campaign calls once, with the SAME config
object/closure it then uses for the negative sweep.

Usage
-----
A campaign wraps its solver as a zero-argument closure that returns a result for
a known-positive control, then asserts re-discovery *before* trusting any
negative::

    from cyclerfinder.verify.positive_control import assert_positive_control

    def run_control() -> SolveResult:
        # IDENTICAL pipeline config to the negative sweep below.
        return solve(known_aldrin_entry, **sweep_config)

    assert_positive_control(
        run_control,
        found=lambda res: res.closes,
        label="Aldrin classic k=1 through the #N sweep config",
    )
    # ... only now run the negative sweep with the same sweep_config ...

The point is that ``run_control`` and the negative sweep share ``sweep_config``
by construction, so a positive control that passes certifies the *pipeline*, not
just the math.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

_T = TypeVar("_T")


class PositiveControlError(AssertionError):
    """A known-positive control was NOT re-found through the campaign's pipeline.

    Raised by :func:`assert_positive_control`. Subclasses :class:`AssertionError`
    so it reads as a test/gate failure, while remaining catchable distinctly by
    campaign drivers that want to annotate the un-trusted negative.
    """


def assert_positive_control(
    run_control: Callable[[], _T],
    *,
    found: Callable[[_T], bool],
    label: str,
) -> _T:
    """Assert a negative campaign re-finds a KNOWN solution before being trusted.

    Runs ``run_control()`` (which MUST use the identical pipeline configuration
    the subsequent negative sweep uses) and asserts ``found(result)`` is true.
    Returns the control result on success so the caller can log it; raises
    :class:`PositiveControlError` otherwise.

    Parameters
    ----------
    run_control:
        Zero-argument closure that runs the method on a known-positive control
        through the SAME config as the negative sweep. Sharing the config by
        construction is the whole point — it certifies the pipeline, not just the
        underlying math.
    found:
        Predicate returning ``True`` when ``run_control``'s result counts as a
        re-discovery of the control (e.g. ``lambda r: r.closes``).
    label:
        Human-readable description of the control + config, quoted in the error.

    Returns
    -------
    The control result (when re-discovery succeeded).

    Raises
    ------
    PositiveControlError
        If ``found(run_control())`` is false — the method failed to re-find a
        solution it should have, so any negative result from the same pipeline is
        uninformative (the #180-class trap) and must not be trusted.
    """
    result = run_control()
    if not found(result):
        raise PositiveControlError(
            f"positive control NOT re-found through the identical pipeline config "
            f"({label}); the method is mis-configured, so any negative result from "
            f"this pipeline is uninformative and must not be trusted (#180 class). "
            f"Fix the config until the control is rediscovered before trusting the "
            f"negative."
        )
    return result


__all__ = ["PositiveControlError", "assert_positive_control"]
