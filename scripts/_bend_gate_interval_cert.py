"""#625 -- shared interval-arithmetic bend-gate certification helper.

Generalizes ``scripts/certify_610_proteus_bend_interval.py``'s Proteus-only
proof-of-concept into a small, reusable, parameterized building block, so
every other ``data/empty_regions.jsonl`` entry that fails its negative for
the SAME reason (the ``#324`` physical max-bend gate, undersized-moon-GM
branch -- :mod:`cyclerfinder.search.physical_sanity`) can be certified with
one function call instead of a bespoke script.

This module carries ONLY the reusable primitives (rigorous interval arcsin,
the interval enclosure of the Bate-Mueller-White max-bend formula, and a
thin structured wrapper around it). It intentionally does NOT itself target
any specific entry -- see ``scripts/certify_625_bend_gate_registry.py`` for
the per-entry driver, and ``scripts/certify_610_proteus_bend_interval.py``
(which now imports :func:`rigorous_arcsin` / :func:`bend_deg_interval` from
here instead of defining them locally -- byte-for-byte the same math, just
de-duplicated) for the original single-entry proof of concept.

Method / library -- see ``certify_610_proteus_bend_interval.py``'s own module
docstring for the full citation trail (Moore/Kearfott/Cloud interval
exclusion tests; ``mpmath.iv`` directed-rounding interval context; the
``atan2``/``sqrt`` rigorous-arcsin identity). Not repeated here verbatim to
avoid drift between two copies of the same prose -- that docstring is the
canonical description of the METHOD; this module is the canonical location
of the CODE.

Honest scope limit (inherited from #610, restated per #625's own task spec):
:func:`certify_bend_gate_over_box` certifies ONLY that no point in the given
``(r_p, V_inf)`` box can clear the physical-bend sub-gate -- it says nothing
about whether the residual/Lambert-closure half of any search that box is
drawn from is exhaustive. Bounding a multi-revolution, branch-selecting
universal-variable Lambert solve rigorously remains a genuinely unresolved
technical obstacle for this project, not something this module attempts.
Every caller MUST carry this caveat forward into whatever it records.
"""

from __future__ import annotations

from typing import Any

try:
    import mpmath as mp  # noqa: F401 -- imported only to probe availability (HAVE_MPMATH)

    HAVE_MPMATH = True
except ImportError:  # pragma: no cover - exercised by the skip-clean test
    HAVE_MPMATH = False


def rigorous_arcsin(iv: Any, t: Any) -> Any:
    """Rigorous interval arcsin via ``atan2(t, sqrt(1 - t**2))``.

    Valid for ``t`` (an ``mpmath.iv`` interval) contained in ``[-1, 1]``.
    ``mpmath.iv`` does not expose ``asin``/``acos`` directly, but does expose
    ``sqrt`` and ``atan2`` as directed-rounding interval primitives -- this
    is the standard identity connecting them (both sides equal, for
    ``t in [0, 1]``, the unique angle in ``[0, pi/2]`` with that sine).
    """
    return iv.atan2(t, iv.sqrt(1 - t**2))


def bend_deg_interval(
    iv: Any,
    mu_km3_s2: float,
    rp_lo_km: float,
    rp_hi_km: float,
    vinf_lo_kms: float,
    vinf_hi_kms: float,
) -> Any:
    """Rigorous interval enclosure of ``max_bend`` (degrees) over a box.

    Box: ``r_p in [rp_lo_km, rp_hi_km]`` km, ``V_inf in [vinf_lo_kms,
    vinf_hi_kms]`` km/s. ``mu_km3_s2`` is treated as an exact point (the
    registry's stored double). Mirrors
    :func:`cyclerfinder.core.flyby.max_bend`'s formula exactly:
    ``sin(delta/2) = mu / (mu + rp * vinf**2)``, ``delta = 2*arcsin(...)``.

    ``max_bend`` is monotonically DECREASING in both ``r_p`` and ``V_inf``
    (larger periapsis or faster flyby both shrink the deflection), so the
    enclosure's supremum is realized exactly at the box's lower corner
    ``(rp_lo_km, vinf_lo_kms)`` -- the upper bounds ``rp_hi_km``/
    ``vinf_hi_kms`` only need to be *some* value at or above the lower bound
    for a mathematically valid box; they do not affect ``sup(bend)``. This
    lets a caller safely widen the upper bounds for a more conservative-
    looking box without changing the certified result.
    """
    mu = iv.mpf(mu_km3_s2)
    rp = iv.mpf([rp_lo_km, rp_hi_km])
    vinf = iv.mpf([vinf_lo_kms, vinf_hi_kms])
    ratio = mu / (mu + rp * vinf**2)
    half_delta = rigorous_arcsin(iv, ratio)
    delta_rad = 2 * half_delta
    return delta_rad * (180 / iv.pi)


def certify_bend_gate_over_box(
    iv: Any,
    *,
    gm_moon_km3_s2: float,
    rp_lo_km: float,
    rp_hi_km: float,
    vinf_lo_kms: float,
    vinf_hi_kms: float,
    gate_deg: float,
    label: str = "",
) -> dict[str, Any]:
    """Certify (or fail to certify) that a bend gate is unsatisfiable over a box.

    Thin, structured wrapper around :func:`bend_deg_interval`: computes the
    rigorous enclosure, extracts its supremum, and compares against
    ``gate_deg``. Returns a dict suitable for embedding directly as a new
    ``empty_regions.jsonl`` field (task #625) -- NOT a full
    :class:`~cyclerfinder.data.empty_regions.EmptyRegionReport`, just this
    one sub-claim.

    Honest scope note (mirrored into the returned dict's own ``scope_note``
    so it travels with the data, not just this module's docstring): this
    certifies ONLY that no ``(r_p, V_inf)`` in the given box can clear the
    physical max-bend sub-gate. It says nothing about whether the box itself
    covers every point a real search's Lambert/residual closure half could
    reach -- that half stays uncertified (see module docstring).
    """
    enclosure = bend_deg_interval(iv, gm_moon_km3_s2, rp_lo_km, rp_hi_km, vinf_lo_kms, vinf_hi_kms)
    sup_bend_deg = float(enclosure.b)
    inf_bend_deg = float(enclosure.a)
    certified = sup_bend_deg < gate_deg
    return {
        "label": label,
        "box": {
            "gm_km3_s2": gm_moon_km3_s2,
            "rp_km": [rp_lo_km, rp_hi_km],
            "vinf_kms": [vinf_lo_kms, vinf_hi_kms],
        },
        "inf_bend_deg": inf_bend_deg,
        "sup_bend_deg": sup_bend_deg,
        "gate_deg": gate_deg,
        "certified": certified,
        "method": "mpmath.iv interval arithmetic (rigorous_arcsin via atan2/sqrt identity)",
        "scope_note": (
            "Certifies only that no (r_p, V_inf) in the recorded box can clear the "
            "physical max-bend sub-gate (the #324 gate). Does NOT certify the "
            "residual/Lambert-closure half of the underlying search is exhaustive; "
            "bounding a multi-rev branch-selecting Lambert solve rigorously remains "
            "unsolved. See task #610/#625."
        ),
    }


__all__ = [
    "HAVE_MPMATH",
    "bend_deg_interval",
    "certify_bend_gate_over_box",
    "rigorous_arcsin",
]
