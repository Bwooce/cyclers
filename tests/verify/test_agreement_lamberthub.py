"""Axis-A lamberthub-path integration test — the Forge, Phase 2 (slow).

Path (a) of :func:`cyclerfinder.verify.agreement.crosscheck_code_paths`
reuses the M7 spec §14 V1 cross-check
(:func:`cyclerfinder.verify.crosscheck.crosscheck_cycler`): every leg
re-solved with ``lamberthub``'s ``izzo2015`` + ``gooding1990`` and
compared to the in-house Lambert. Marked slow because it imports the
astropy ephemeris and runs three Lambert solvers per leg.
"""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler
from cyclerfinder.search.phase_match import phase_signature_from_catalogue_entry
from cyclerfinder.verify.agreement import crosscheck_code_paths
from cyclerfinder.verify.crosscheck import V1_TOLERANCE_MPS
from cyclerfinder.verify.real_closure import (
    _resolve_real_t_start,
    construct_real_ephemeris_cycler,
)
from tests.data._catalogue_loader_m6b import load_m6b_entries

ALDRIN_PRIORITY = datetime(1985, 10, 28, tzinfo=UTC)


@pytest.fixture(scope="module")
def astropy_ephem() -> Ephemeris:
    return Ephemeris(model="astropy")


@pytest.fixture(scope="module")
def aldrin_real_cycler(astropy_ephem: Ephemeris) -> Cycler:
    """Real-ephemeris Aldrin cycler (single-rev E->M->E Lambert chain)."""
    entry = next(e for e in load_m6b_entries() if e["id"] == "aldrin-classic-em-k1-outbound")
    signature = phase_signature_from_catalogue_entry(entry)
    t_start = _resolve_real_t_start(signature, astropy_ephem, ALDRIN_PRIORITY)
    assert t_start is not None
    return construct_real_ephemeris_cycler(entry, astropy_ephem, t_start)


@pytest.fixture(scope="module")
def aldrin_inbound_real_cycler(astropy_ephem: Ephemeris) -> Cycler:
    """Real-ephemeris Aldrin INBOUND cycler (single-rev M->E->M Lambert chain).

    Built like-for-like with the outbound fixture above: same loader, same
    phase-signature resolution, same constructor — only the catalogue row
    differs (the inbound / down-escalator twin). This is the row Part 1 of
    the validation campaign (#125) promotes to spec §14 V1.
    """
    entry = next(e for e in load_m6b_entries() if e["id"] == "aldrin-classic-em-k1-inbound")
    signature = phase_signature_from_catalogue_entry(entry)
    t_start = _resolve_real_t_start(signature, astropy_ephem, ALDRIN_PRIORITY)
    assert t_start is not None
    return construct_real_ephemeris_cycler(entry, astropy_ephem, t_start)


@pytest.mark.slow
def test_report_includes_lamberthub_path(
    aldrin_real_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """The report includes the reused lamberthub cross-check and it passes.

    Integration check that path (a) is wired into the combiner and that,
    with an ephemeris supplied, the lamberthub path becomes available and
    every leg agrees to < :data:`V1_TOLERANCE_MPS`.
    """
    report = crosscheck_code_paths(aldrin_real_cycler, astropy_ephem)
    path_a = report.lamberthub
    assert path_a.available is True
    assert path_a.per_leg, "no legs cross-checked"
    assert path_a.passed is True, (
        f"lamberthub disagreement {path_a.max_diff_mps} m/s >= {V1_TOLERANCE_MPS} m/s"
    )
    assert path_a.max_diff_mps < V1_TOLERANCE_MPS


@pytest.mark.slow
def test_real_eph_paths_a_and_c_pass_b_flags_model_mismatch(
    aldrin_real_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """All three Axis-A paths run on the real-ephemeris Aldrin cycler;
    the lamberthub (a) and Kepler re-propagation (c) paths PASS, while the
    circular-coplanar resonance construction (b) correctly FLAGS the
    model mismatch.

    This is the honest real-ephemeris behaviour and a useful teeth check
    on the combiner's veto rule. Path (b) is a *circular-coplanar*
    single-ellipse construction; the real-ephemeris Aldrin geometry is a
    3-encounter E-M-E chain whose per-encounter V_inf (~14.5 / 13.2 km/s)
    departs sharply from the coplanar single-ellipse construction
    (~6.4 / 8.3 km/s) because real Mars eccentricity and the multi-
    encounter topology break the single-ellipse assumption. The path is
    therefore *available but failing*, which the combiner treats as a veto:
    ``agreed`` is False even though two independent paths (a, c) agree.
    The clean construction-vs-construction agreement is exercised on the
    circular-coplanar seed in ``test_agreement.py``.
    """
    report = crosscheck_code_paths(aldrin_real_cycler, astropy_ephem)
    assert report.lamberthub.available is True
    assert report.lamberthub.passed is True
    assert report.kepler_reprop.available is True
    assert report.kepler_reprop.passed is True
    # All three paths run; path (b) flags the coplanar-vs-real mismatch.
    assert report.n_paths_available == 3
    assert report.construction_optimiser.available is True
    assert report.construction_optimiser.passed is False
    assert report.n_paths_passed == 2
    # A failing available path vetoes the overall verdict.
    assert report.agreed is False


@pytest.mark.slow
def test_inbound_real_eph_lamberthub_and_kepler_paths_pass(
    aldrin_inbound_real_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """Part 1 (#125): the real-DE440 Aldrin INBOUND cycler clears spec §14 V1.

    Same evidence the outbound row earns its V1 from — the two §14 V1 halves:
    every leg re-solved with ``lamberthub`` izzo2015 + gooding1990 agrees to
    < :data:`V1_TOLERANCE_MPS` (path a), AND the Kepler forward
    re-propagation residual passes (path c). Both paths are available and
    pass at the same bars as outbound, so the overall verdict ``agreed`` is
    True.

    The honest difference from the outbound twin: the circular-coplanar
    resonance construction (path b) is *unavailable* here, not
    available-but-failing. The inbound chain starts on the short ~146 d
    Mars->Earth arc, whose real-ephemeris Lambert solution yields a
    heliocentric ``(a, e)`` that does not map onto the analytic single-
    ellipse crossing construction, so :func:`_resonant_vinf_for_cycler`
    returns ``None`` (recorded as unavailable, not a failure). The V1
    verdict is unaffected: §14 V1 is defined by the lamberthub-agreement
    and Kepler-re-propagation halves, both of which pass.
    """
    report = crosscheck_code_paths(aldrin_inbound_real_cycler, astropy_ephem)
    # Path (a): lamberthub leg agreement < V1_TOLERANCE_MPS.
    assert report.lamberthub.available is True
    assert report.lamberthub.per_leg, "no legs cross-checked"
    assert report.lamberthub.passed is True, (
        f"lamberthub disagreement {report.lamberthub.max_diff_mps} m/s >= {V1_TOLERANCE_MPS} m/s"
    )
    assert report.lamberthub.max_diff_mps < V1_TOLERANCE_MPS
    # Path (c): Kepler forward re-propagation residual passes.
    assert report.kepler_reprop.available is True
    assert report.kepler_reprop.passed is True
    # Path (b) is unavailable (single-ellipse crossing undefined for the
    # short M->E first leg) — recorded as unavailable, not failing.
    assert report.construction_optimiser.available is False
    # Two independent paths available, both pass -> agreed verdict.
    assert report.n_paths_available == 2
    assert report.n_paths_passed == 2
    assert report.agreed is True
