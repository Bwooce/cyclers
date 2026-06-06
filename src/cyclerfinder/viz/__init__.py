# AP_FLAKE8_CLEAN
"""Optional-extra visualisation package for cyclerfinder.

matplotlib is the optional ``[viz]`` extra — it is **lazy-imported** inside
:mod:`cyclerfinder.viz.plots` (Agg backend) and never a hard dependency. Code
that needs a plot but runs on a base install gets a clear, actionable
:class:`MissingVizExtra` error rather than an opaque ``ImportError``.
"""

from __future__ import annotations


class MissingVizExtra(RuntimeError):  # noqa: N818  (plan-mandated name; not "...Error")
    """Raised when a plotting call needs matplotlib but the [viz] extra is absent.

    Carries the install hint so the message is actionable.
    """


# Import the submodule so ``from cyclerfinder.viz import plots`` works without a
# circular ``from ... import plots`` re-export. plots.py does NOT import
# matplotlib at module level, so this stays cheap on a base install.
from cyclerfinder.viz import plots as plots  # noqa: E402

__all__ = ["MissingVizExtra", "plots"]
