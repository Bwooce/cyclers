"""N-body Phase C.0: the shooter honours the #135 seeding verdict (plan Phase C).

#135 verdict (russell12-likeforlike): the basin problem is SEEDING, not solver.
The shooter must be seeded from the #133 near-miss survey, never a blind scan.
This guard pins that intent in the module docstring so it cannot be lost.
"""

from __future__ import annotations

import inspect


def test_shooter_documents_nearmiss_seeding() -> None:
    from cyclerfinder.nbody import shooter

    doc = inspect.getdoc(shooter) or ""
    assert "near-miss" in doc.lower()
    assert "135" in doc  # cross-reference the verdict that mandated it
