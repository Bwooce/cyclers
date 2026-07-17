"""#610 -- certified non-existence via interval arithmetic (proof-of-concept).

Background (``data/OUTSTANDING.md`` #605/#610). Every entry in
``data/empty_regions.jsonl`` is, by the registry's own stated epistemic
status, a *numerical* negative: a discrete grid/enumeration found no
gate-passing candidate, which is evidence of absence conditional on the
method (grid resolution, corrector, tolerance), never proof of absence. This
script attempts a genuinely different, stronger class of result over ONE
existing entry: a compact-region non-existence CERTIFICATE built from
rigorous interval arithmetic, following the classical interval-analysis
exclusion-test method (Moore, R.E., Kearfott, R.B. & Cloud, M.J.,
*Introduction to Interval Analysis*, SIAM, 2009 -- see esp. Ch. 1-4 on
interval extensions of real functions and range enclosure/exclusion tests).
No interval-Newton root ITERATION is needed here (there is no equation to
solve) -- :func:`max_bend` is a monotonic scalar function of two variables,
so a single rigorous *range enclosure* over a box is already a valid
exclusion test: if the enclosure's supremum is below the gate threshold, NO
point in the (uncountable) box -- not just the discrete grid points actually
evaluated -- can pass the gate. That is the precise sense in which this
upgrades the entry from "conditional on method" to "certified over this
box."

Target entry
------------
``neptune-triton-proteus-symmetric-closure-599-2026-07-15`` (task #599).
Chosen over the other #563-method entries (Titan-Iapetus #575, Galilean
#576/#577, Uranian 3-moon #600) because its own recorded "verdict"/
"interpretation" already isolates a SINGLE, low-dimensional, non-transcendental
sub-condition as the entire reason for the negative: every one of the 104
residual-gate survivors (of 1024 evaluated) fails specifically because
Proteus's GM is too small to deliver >= 5 deg of ballistic bend at the
V_inf regime this symmetric-closure family produces
(:mod:`cyclerfinder.search.physical_sanity`, the #324 gate). That specific
claim -- "no achievable (r_p, V_inf) at Proteus reaches 5 deg of bend" -- is
EXACTLY the closed-form Bate-Mueller-White patched-conic formula

.. math::

    \\sin(\\delta_\\max/2) = \\mu / (\\mu + r_p V_\\infty^2)

(:func:`cyclerfinder.core.flyby.max_bend`), a composition of +, *, /, and
arcsin -- tractable to bound rigorously with interval arithmetic. The OTHER
half of the entry's search (which (rel_offset, tof, n_rev) combinations reach
residual < 0.05 km/s at all) is a multi-revolution Lambert/Kepler-equation
solve -- transcendental, iterative, branch-selecting -- and is NOT attempted
here; see the "What this does NOT certify" section below. That scope split
is the honest, bounded deliverable this task asks for.

Method / library
-----------------
``mpmath.iv`` -- mpmath's built-in interval-arithmetic context (Johansson,
F. et al., *mpmath: a Python library for arbitrary-precision floating-point
arithmetic*, http://mpmath.org, interval context documented at
http://mpmath.org/doc/current/contexts.html#the-interval-context). Verified
directly (not just cited) against ``math.asin`` at 5 sample points before
use (see ``tests/scripts/test_certify_610_proteus_bend_interval.py``).
``mpmath.iv`` does not expose ``asin``/``acos`` directly, so this module
builds a rigorous interval arcsin from two primitives ``iv`` DOES provide
(``sqrt``, ``atan2``) via the standard identity, valid for
``t in [-1, 1]``::

    arcsin(t) = atan2(t, sqrt(1 - t**2))

Honesty note on assurance level: ``mpmath.iv`` is a widely-used
directed-rounding software implementation, not a formally verified one (cf.
INTLAB for MATLAB, or a hardware-directed-rounding-based C library). The
enclosures below were spot-checked against non-interval ``math.asin``/
``core.flyby.max_bend`` at multiple points and the enclosure width is
inspected to confirm it stays microscopic (see test file) -- but this is
still "a careful, checked software computation," not a machine-checked
formal proof. Reported as such.

Result (see module-level docstring in the run's stdout / this file's tests
for the numbers): CERTIFIED -- ``sup(bend_deg)`` stays strictly below the 5
deg gate over BOTH a tight, data-grounded box and a substantially widened
one; see ``main()`` for the exact boxes and numbers.

What this DOES certify
-----------------------
For the Proteus encounter of the #599 entry's symmetric-closure family: NO
periapsis radius in ``[r_safe, r_safe_upper]`` km and NO hyperbolic excess
speed in the certified V_inf interval (an uncountable continuum, not a grid)
can ever deliver >= 5 deg of ballistic bend. This directly and rigorously
grounds the entry's own "Proteus's GM is too small" interpretation as a
continuum-strength claim, not merely "the 104 grid survivors we happened to
enumerate all failed."

What this does NOT certify
----------------------------
1. The residual/closure half of the search (does some (rel_offset, tof,
   n_rev) achieve residual < 0.05 km/s at all) is NOT interval-certified
   here. Bounding a multi-revolution Lambert solve (branch-selecting,
   built on an iterative universal-variable Kepler solve) rigorously is a
   real, unresolved technical obstacle for this proof-of-concept's scope --
   reported honestly per the task's own "genuine technical obstacle is a
   valid result" allowance, not forced.
2. This certifies ONLY the Proteus encounter of ONE entry's ONE failure
   mode. It is not a general framework and makes no claim about Triton's
   bend (which is large and NOT the reason for the entry's negative), nor
   about any other empty_regions.jsonl entry.
3. Registry GM/radius/safe-altitude values are treated as exact (their own
   OBSERVATIONAL uncertainty is a separate concern from the interval-
   arithmetic ROUNDING-rigor claim made here).

No catalogue.yaml / empty_regions.jsonl writes. This is a proof-of-concept
script only (task #610); registering any upgraded claim in the actual
registry is a follow-up decision for the coordinating session.

Update (2026-07-17, task #625): the two reusable primitives this script
built (:func:`rigorous_arcsin`, :func:`bend_deg_interval`) have been
factored out into ``scripts/_bend_gate_interval_cert.py`` and are now
IMPORTED from there rather than defined locally -- byte-for-byte the same
math, de-duplicated so #625's other target entries (and any future one)
reuse the exact same, already-tested implementation instead of a second
hand-copy. This script (and its own test file) is unchanged in behavior;
``certify610.rigorous_arcsin`` / ``certify610.bend_deg_interval`` still
resolve as module attributes (imported names), so existing tests need no
changes. See ``scripts/certify_625_bend_gate_registry.py`` for the
generalized driver that certifies the other entries.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
# Also on sys.path so `scripts.X` (dotted) resolves standalone, matching how
# tests/scripts/*.py already address these modules -- needed because a couple
# of sibling scripts (this one included) are reachable BOTH ways depending on
# caller, and mypy needs a single, consistent module identity to avoid a
# "Source file found twice under different module names" error.
sys.path.insert(0, str(ROOT))

from cyclerfinder.core.satellites import SATELLITES  # noqa: E402
from cyclerfinder.search.physical_sanity import DEFAULT_MIN_USEFUL_BEND_DEG  # noqa: E402

try:
    import mpmath as mp

    HAVE_MPMATH = True
except ImportError:  # pragma: no cover - exercised by the skip-clean test
    HAVE_MPMATH = False

# #625: the two reusable primitives now live in a shared module -- imported
# (not redefined) so this script and #625's driver share one implementation.
# Re-exported under their original names so `certify610.rigorous_arcsin` /
# `certify610.bend_deg_interval` keep resolving exactly as before (existing
# test file imports them off this module, unchanged). `__all__` makes the
# re-export explicit for mypy's strict-mode implicit-reexport check.
from scripts._bend_gate_interval_cert import bend_deg_interval, rigorous_arcsin  # noqa: E402

__all__ = ["bend_deg_interval", "rigorous_arcsin"]


def _proteus_subgate_vinf_range() -> tuple[float, float, int, int]:
    """Reproduce the real #599 Proteus V_inf range among residual-gate survivors.

    Reuses (does not modify) the #558/#563/#599 machinery verbatim, exactly
    as the enumeration script itself does, to ground the certified box in
    the entry's OWN data rather than a hand-copied prose number.

    Returns ``(vinf_min_kms, vinf_max_kms, n_evaluated, n_subgate)``.
    """
    import itertools

    from scan_558_uranus_all_pairs_offset_sweep import (
        GATE_RESIDUAL_KMS,
        gate_candidate,
        residual_at_point,
    )

    from scripts.enumerate_563_symmetric_closures import N_REV_VALUES, REL_OFFSETS_DEG, pair_n_max

    primary = "Neptune"
    directions = [("Triton", "Proteus"), ("Proteus", "Triton")]
    vinfs_at_proteus: list[float] = []
    n_total = 0
    n_subgate = 0
    for anchor, flyby in directions:
        t_syn, p_a, p_b, n_max = pair_n_max(anchor, flyby, primary=primary)
        sqrt_papb = (p_a * p_b) ** 0.5
        for n in range(1, n_max + 1):
            target_tof_days = n * t_syn / 2.0
            target_tof_scale = target_tof_days / sqrt_papb
            for n0, n1 in itertools.product(N_REV_VALUES, N_REV_VALUES):
                for rel in REL_OFFSETS_DEG:
                    n_total += 1
                    pt = residual_at_point(
                        anchor,
                        flyby,
                        rel_offset_deg=rel,
                        tof_scale=target_tof_scale,
                        n_rev=(n0, n1),
                        primary=primary,
                    )
                    if pt is None or pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                        continue
                    n_subgate += 1
                    gated = gate_candidate(anchor, flyby, pt, primary=primary)
                    vinfs = gated["vinf_per_encounter_kms"]  # [anchor, flyby, anchor]
                    proteus_idx = 1 if flyby == "Proteus" else 0
                    vinfs_at_proteus.append(vinfs[proteus_idx])

    return min(vinfs_at_proteus), max(vinfs_at_proteus), n_total, n_subgate


def main() -> int:
    if not HAVE_MPMATH:
        print(
            "[610] mpmath not installed -- run `uv run --with mpmath python "
            "scripts/certify_610_proteus_bend_interval.py`, or `uv sync --extra interval`.",
            flush=True,
        )
        return 1

    mp.mp.dps = 50
    iv = mp.iv
    iv.dps = 50

    proteus = SATELLITES["Proteus"]
    mu = proteus.mu_km3_s2
    rp_safe = proteus.radius_eq_km + proteus.safe_alt_km
    gate_deg = DEFAULT_MIN_USEFUL_BEND_DEG

    print(
        f"[610] Proteus: mu={mu} km^3/s^2, r_p_safe={rp_safe} km, gate={gate_deg} deg", flush=True
    )

    vinf_min, vinf_max, n_total, n_subgate = _proteus_subgate_vinf_range()
    print(
        f"[610] reproduced #599 data: {n_total} evaluated, {n_subgate} pass residual "
        f"sub-gate; Proteus V_inf range among survivors = [{vinf_min:.6f}, {vinf_max:.6f}] km/s",
        flush=True,
    )

    # Box A: tight, data-grounded -- exactly the V_inf range that occurs among
    # every one of #599's 104 residual-gate survivors, r_p pinned at the
    # registry safety floor (the periapsis choice that MAXIMIZES bend; any
    # larger r_p only decreases it further, per the formula's monotonicity).
    bend_a = bend_deg_interval(iv, mu, rp_safe, rp_safe, vinf_min, vinf_max)
    sup_a = float(bend_a.b)
    print(
        f"[610] Box A (data-grounded)  r_p=[{rp_safe},{rp_safe}] km, "
        f"V_inf=[{vinf_min:.6f},{vinf_max:.6f}] km/s -> bend_deg enclosure = {bend_a}, "
        f"sup={sup_a:.6f}",
        flush=True,
    )
    certified_a = sup_a < gate_deg
    print(f"[610] Box A certified (sup < {gate_deg} deg): {certified_a}", flush=True)

    # Box B: widened, conservative -- r_p over a full order of magnitude above
    # the safety floor, V_inf down to a rigorously-derived critical margin
    # (analytically, bend=5deg at V_inf~0.429 km/s for r_p=rp_safe; 0.45 km/s
    # is a safety margin below the true 104-survivor minimum of ~1.82 km/s by
    # a factor of >4x) up through a generous 20 km/s ceiling.
    rp_hi = 10.0 * rp_safe
    vinf_lo_margin = 0.45
    vinf_hi_wide = 20.0
    bend_b = bend_deg_interval(iv, mu, rp_safe, rp_hi, vinf_lo_margin, vinf_hi_wide)
    sup_b = float(bend_b.b)
    print(
        f"[610] Box B (widened)       r_p=[{rp_safe},{rp_hi}] km, "
        f"V_inf=[{vinf_lo_margin},{vinf_hi_wide}] km/s -> bend_deg enclosure = {bend_b}, "
        f"sup={sup_b:.6f}",
        flush=True,
    )
    certified_b = sup_b < gate_deg
    print(f"[610] Box B certified (sup < {gate_deg} deg): {certified_b}", flush=True)

    if certified_a and certified_b:
        print(
            "[610] RESULT: CERTIFIED non-existence of any Proteus flyby bend >= "
            f"{gate_deg} deg over both boxes above -- upgrades the #599 entry's "
            "'Proteus GM too small' interpretation from a 104-point grid finding "
            "to a continuum-strength interval-arithmetic certificate over this "
            "sub-region of the entry (the physical-bend sub-gate only; the "
            "residual/Lambert-closure half of the search remains uncertified, "
            "see module docstring).",
            flush=True,
        )
        return 0
    print("[610] RESULT: NOT certified over one or both boxes -- see enclosures above.", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
