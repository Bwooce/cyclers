"""INDEPENDENT n-body cross-check of the reachable App-C V3-batch parents (#170).

Scales the proven S1L1 gate (``test_s1l1_corrected_nbody.py``) to the two
``russell-ch4-*`` V0 parents that carry a sourced App-C "DATA NECESSARY TO REPRODUCE"
block (Russell 2004 Appendix C): ``8.049gGf2`` (#188) and ``8.165Gfh-f2`` (#192).
The other seven V0 ``russell-ch4-*`` rows are NOT-REACHABLE (no sourced App-C parent
number / block) and are recorded+skipped in the batch note, not forced.

VERDICT (pinned): **PARTIAL** for both. The Mars encounters reconstruct cleanly and
land in-band on an INDEPENDENT integrator at the published per-leg v_inf — but neither
is V3-grade: both are POWERED cyclers whose OWN published App-C total Δv (420 / 1678
m/s over 7 cycles, Table 5.5 = the ``block.total_dv_kms`` field) is far over the spec
§14 V3 maintenance budget (120 m/s). That sourced total is the BINDING V3-rejection
here — it is Russell's number, epoch- and model-independent. So: NO V3 writeback. The
encounter geometry is real, the fuel cost is not V3-grade.

The Sun-only patched-conic continuous-from-one-seed TCM (the #169 proxy) is recorded
as corroborating EVIDENCE only: it is a LOWER bound on the real maintenance and is NOT
the gate. After the #198 63 s UTC/TDB epoch fix (commit 439d279) it re-derives to
114.4 m/s (#188) and 2020.7 m/s (#192) — down from the pre-fix 163.6 / 2040.6 m/s
recorded in ``docs/notes/2026-06-08-appc-v3-batch-results.md``. Note #188's proxy now
sits UNDER the 120 m/s budget (114.4) even though the published powered total is 420
m/s — exactly why the verdict is grounded on the sourced total, not the Sun-only proxy.
A faithful continuous measure needs the Mars-perturbed run with B-plane flyby targeting
(M7); that is not this test's job. See
``docs/notes/2026-06-23-appc-s1l1-tcm-epoch-rederivation.md``.

What independence means (same as the S1L1 gate): the only thing seeded is each leg's
departure state (App-C v_inf + the real DE440 planet velocity at the printed epoch +
the printed mid-leg Δv). Where the spacecraft goes / where DE440 Mars is / the
achieved miss+v_inf is pure integrator+ephemeris output, never fit.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.constants import (  # noqa: E402
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
)
from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.search.appc_corrected import (  # noqa: E402
    REACHABLE_BLOCKS,
    AppCArc,
    AppCBlock,
    build_seeded_arcs,
    continuous_chain,
)

# Same 3-Mars-SOI confirmation band as S1L1 / #165 — kept IDENTICAL, never loosened.
_MARS_SOI_AU = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
_ENCOUNTER_MISS_TOL_AU = 3.0 * _MARS_SOI_AU  # ~0.0116 AU
_V3_TCM_BUDGET_KMS = 0.120  # spec §14 V3 maintenance budget

# Computed regression anchors (NOT validation goldens, NOT the verdict gate): the
# Sun-only patched-conic continuous-from-one-seed TCM (#169 proxy) under the
# CORRECTED ephemeris (post-#198 63 s UTC/TDB epoch fix, commit 439d279). These pin
# "the code still computes the same proxy"; the V3-rejection is grounded on the
# sourced ``block.total_dv_kms`` (420 / 1678 m/s). The proxy is a lower bound — for
# #188 it sits UNDER 120 m/s, which is the whole reason it is not the gate. See
# docs/notes/2026-06-23-appc-s1l1-tcm-epoch-rederivation.md.
_SUNONLY_PROXY_TCM_MPS: dict[str, float] = {
    "russell-ch4-8.049gGf2": 114.4,
    "russell-ch4-8.165Gfh-f2": 2020.7,
}


def _propagate_with_dv_nbody(
    arc: AppCArc, ephem: Ephemeris, prop: RestrictedNBody
) -> tuple[float, float, bool]:
    """Sun-only IAS15 propagation of one Mars-transit leg to arrival, applying mid-leg Δv.

    Returns ``(miss_au, vinf_kms, converged)`` against real DE440 Mars. Sun-only is
    Russell's patched-conic cruise model (IAS15 pinned to analytic two-body by GOLDEN
    GATE 1), the apples-to-apples independent reproduction of the published seed.
    """
    converged = True
    if arc.t_dv_sec is not None and arc.dv_km_s is not None and arc.t_dv_sec < arc.t1_sec:
        leg1 = prop.propagate(
            arc.r0_km, arc.v0_km_s, arc.t0_sec, arc.t_dv_sec, accuracy=1e-11, max_wall_sec=120.0
        )
        converged = converged and bool(leg1.converged)
        r_dv = np.asarray(leg1.r_km, dtype=np.float64)
        v_dv = np.asarray(leg1.v_km_s, dtype=np.float64) + arc.dv_km_s
        out = prop.propagate(
            r_dv, v_dv, arc.t_dv_sec, arc.t1_sec, accuracy=1e-11, max_wall_sec=120.0
        )
    else:
        out = prop.propagate(
            arc.r0_km, arc.v0_km_s, arc.t0_sec, arc.t1_sec, accuracy=1e-11, max_wall_sec=120.0
        )
    converged = converged and bool(out.converged)
    r_m, v_m = ephem.state("M", arc.t1_sec)
    miss_au = float(np.linalg.norm(out.r_km - np.asarray(r_m)) / AU_KM)
    vinf_kms = float(np.linalg.norm(out.v_km_s - np.asarray(v_m)))
    return miss_au, vinf_kms, converged


@pytest.mark.slow
@pytest.mark.parametrize("catalogue_id", sorted(REACHABLE_BLOCKS))
def test_appc_batch_encounters_in_band_but_published_dv_over_v3_budget(
    catalogue_id: str,
) -> None:
    """Independent n-body: App-C parents reconstruct in-band; published Δv > V3 budget.

    PINNED VERDICT: PARTIAL. (1) all 7 Mars encounters land in-band on an INDEPENDENT
    Sun-only IAS15 at the published per-leg v_inf (the encounter geometry IS real);
    (2) the SOURCED published App-C total Δv (``block.total_dv_kms`` = 420 / 1678 m/s,
    Table 5.5) exceeds the 120 m/s V3 budget — these are powered cyclers. That sourced
    total is the binding V3-rejection (epoch- and model-independent). Therefore NO V3
    writeback. The Sun-only continuous-from-one-seed TCM (#169 proxy) is checked only
    as corroborating EVIDENCE / a regression anchor — it is a lower bound and for #188
    now reads UNDER 120 m/s after the #198 epoch fix, which is exactly why it is not
    the gate. If a future change moves an encounter out of band, drops the published
    total under budget, or drifts the proxy off its corrected anchor, this asserts
    deliberately — re-derive, do not loosen the band or the budget."""
    block: AppCBlock = REACHABLE_BLOCKS[catalogue_id]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # astropy ERFA dubious-year at far-future epochs
        ephem = Ephemeris("astropy")
        prop = RestrictedNBody("rebound")
        arcs = [a for a in build_seeded_arcs(block, ephem) if a.is_mars_transit]
        assert len(arcs) == 7, f"{catalogue_id}: expected 7 Mars-transit legs, got {len(arcs)}"

        # (1) Encounter geometry — independent integrator, in-band at App-C v_inf.
        for arc in arcs:
            arrival_no = arc.leg_no + 1
            _pub_tof, pub_vinf = block.mars_transit[arrival_no]
            miss_au, vinf_kms, converged = _propagate_with_dv_nbody(arc, ephem, prop)
            assert converged, f"{catalogue_id} M{arrival_no}: Sun-only leg did not converge"
            assert miss_au < _ENCOUNTER_MISS_TOL_AU, (
                f"{catalogue_id} M{arrival_no}: miss {miss_au:.2e} AU outside the "
                f"{_ENCOUNTER_MISS_TOL_AU:.4f} AU band"
            )
            assert vinf_kms == pytest.approx(pub_vinf, abs=5e-3), (
                f"{catalogue_id} M{arrival_no}: v_inf {vinf_kms:.4f} vs published {pub_vinf}"
            )

        # (2) Sun-only patched-conic continuous-from-one-seed TCM (#169 proxy).
        nodes = continuous_chain(block, ephem)
        tcm_kms = sum(n.dv_total_kms for n in nodes)

    # BINDING V3-rejection: the SOURCED published App-C total Δv (Russell Table 5.5)
    # is far over the 120 m/s budget — these are powered cyclers. Epoch-independent.
    assert block.total_dv_kms > _V3_TCM_BUDGET_KMS, (
        f"{catalogue_id}: published App-C total Δv {block.total_dv_kms * 1000:.0f} m/s "
        f"should exceed the 120 m/s V3 budget (these are powered parents)."
    )

    # CORROBORATING EVIDENCE (regression anchor, NOT the gate): the Sun-only proxy is a
    # lower bound on the real maintenance. Pinned to its corrected-ephemeris value; for
    # #188 this is UNDER 120 m/s, which is precisely why the gate above uses the sourced
    # published total, not this proxy. Drift off the anchor => re-derive (see note).
    assert tcm_kms * 1000.0 == pytest.approx(_SUNONLY_PROXY_TCM_MPS[catalogue_id], abs=1.0), (
        f"{catalogue_id}: Sun-only proxy TCM {tcm_kms * 1000:.1f} m/s drifted from the "
        f"corrected-epoch anchor {_SUNONLY_PROXY_TCM_MPS[catalogue_id]} m/s "
        f"(post-#198 epoch fix). Re-derive; do not silently re-pin."
    )
