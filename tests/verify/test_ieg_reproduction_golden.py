"""#480 M1 reproduction golden — Hernandez-Jones-Jesick 2017 IEG triple cycler (AAS 17-608).

The reproduction VERDICT is a characterized NEGATIVE: the EGGIE 4-synodic ballistic
cycler does not reconverge in real Galilean ephemeris from the published-invariant seed
(it relaxes into an off-paper high-energy basin). See
``docs/notes/2026-06-27-480-ieg-reproduction-verdict.md`` for the full numbers, and
``tests/data/test_v4_jupiter.py::test_ieg_jovian_shoot_at_best_epoch`` for the
affirmative, executed not-converged assertion.

So this file does two things:

* keeps the **convergence-independent** citation gate ACTIVE — the golden's EXPECTED side
  is anchored on a verified-against-source citation (#484/#486), guarding the
  concept-collision trap (memory feedback_ground_citations_against_content); and
* SKIPS the reproduction-match assertion, because there is no converged tour to match. The
  sourced tolerances are recorded here unchanged (NOT loosened — "math decides"); the skip
  documents the publication gap rather than forcing a pass.
"""

from __future__ import annotations

import pytest

from cyclerfinder.search.literature_check import anchor_for_key

# SOURCED from docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
# (Hernandez-Jones-Jesick 2017, AAS 17-608, Table 4 — EGGIE 4-synodic ballistic cycler).
# These EXPECTED values trace to the paper, never to code (feedback_golden_tests_sourced_only).
EGGIE_TSYN_DAYS = 7.05
EGGIE_TOTAL_TOF_DAYS = 28.22
EGGIE_GANYMEDE_VINF_KMS = 7.07  # both Ganymede flybys equal
EGGIE_TOTAL_DV_MS = 0.70  # near-ballistic
LAPLACE_RATIO = (1, 2, 4)  # Ganymede:Europa:Io


def test_ieg_citation_is_decision_grade() -> None:
    """The reproduction target resolves to a verified-against-source Jovian IEG anchor.

    Convergence-independent: this guards the citation identity (the #480 hallucination
    trap was attaching the Galilean claim to the VEM sibling AAS 17-577). It must pass
    regardless of whether the tour reproduces.
    """
    a = anchor_for_key("hernandez-2017-ieg-608")  # raises if unresolved/unverified (#484/#486)
    assert a.system == "jovian"
    assert {"Io", "Europa", "Ganymede"} <= set(a.body_set)
    assert a.provenance == "verified-against-source"


@pytest.mark.skip(
    reason="#480 EGGIE not yet reproduced as a ballistic n-body cycler. M1: off-paper basin "
    "(dV ~5.9 km/s). Follow-up 1 SOLVED the basin (resonant-conic seed, all 3 V∞ on Table 4) "
    "and built the corrector program — analytic state+STM co-integrator + block-bidiagonal "
    "Jacobian (~40x faster than FD) + sub-arc multiple shooting — but a ~0.06-0.1 km/s "
    "velocity-continuity wall (localized to the Io perijove) is ROBUST across all four "
    "correctors (FD / analytic-STM / epoch-free / sub-arc): a basin/model feature, not a "
    "numerics one. See docs/notes/2026-06-29-480-eggie-stage4-subarc-verdict.md (+ stage2/3 "
    "+ M1 verdicts). CORRECTION (forward-verify): the resonant-conic seed is structurally "
    "NON-ballistic — the construction pins departure V∞ to the conic but not the equal-in/out "
    "V∞ ballistic-flyby constraint, so Ganymede #2 collapses to ~4 km/s (vs 7.07) and the Io "
    "flyby goes sub-surface; spike-4's analytic-V∞ match was necessary-not-sufficient, and the "
    "3-D-B-plane hypothesis is withdrawn (paper model is coplanar). See "
    "docs/notes/2026-06-29-480-eggie-forward-verify-correction.md. Un-skip only when a true "
    "ballistic-cycler construction (equal-V∞ + bend-feasible + periodic) yields a converged "
    "0.70 m/s EGGIE — WITHOUT loosening the tolerances below."
)
def test_eggie_reproduction_matches_published_invariants() -> None:
    """Reproduction gate (currently un-runnable — no converged tour).

    The sourced tolerances are preserved verbatim for the day a converged seed exists:
    Ganymede V∞ within abs=0.5 km/s of 7.07; correction ΔV < 50 m/s (EIGE-class real-eph
    maintenance ceiling, digest pp.10-11). Do NOT relax these to force a pass.

    Corrector path is ``jovian_shoot`` (Jupiter-central) — Task 4 proved the heliocentric
    ``shooter.shoot()`` raises KeyError on 'Io'; this is the verified Jovian corrector.
    The seed is taken at its best epoch (the multi-rev V∞-matched seed), not the paper
    epoch, since the epoch scan found a sharp seed-defect minimum at -4.54 d.
    """
    from cyclerfinder.nbody.jovian import jovian_shoot
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed

    # Best epoch from scripts/_ieg_rerun_scan_480.py (-4.54 d from the 02-Oct-2020 paper epoch).
    best_departure_et = 654519744.0
    res = jovian_shoot(
        ieg_eggie_seed(departure_et=best_departure_et),
        max_nfev=80,
        max_wall_sec=600.0,
    )
    assert res.converged, f"IEG EGGIE did not converge in real ephemeris: defect={res.defect_norm}"
    gany_vinf = [
        v
        for moon, v in zip(res.sequence, res.vinf_per_encounter_kms, strict=False)
        if moon == "Ganymede"
    ]
    for v in gany_vinf:
        assert v == pytest.approx(EGGIE_GANYMEDE_VINF_KMS, abs=0.5)
    assert res.correction_dv_kms * 1000.0 < 50.0  # m/s
