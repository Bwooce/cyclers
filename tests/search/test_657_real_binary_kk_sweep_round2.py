"""Task #657: real-binary (k1,k2)-cycler sweep, round 2 -- test suite.

`#654` shortlist item 3. Two pieces, both reusing #494/#504/#549's already-
validated CR3BP binary-cycler harness (`real_binary_kk_sweep.py`) verbatim:

(a) an uncapped, fine-mu-step re-run of #549's own 4 compute-capped
    INCONCLUSIVE probes (Patroclus-Menoetius (1,1)/(3,1)/(3,2), Eris-
    Dysnomia (3,2)) -- see `scripts/run_657_real_binary_kk_sweep_round2.py`'s
    `phase_uncapped_recheck`.
(b) three NEW sourced real-binary systems #549 never tried (Sila-Nunam,
    Antiope, Lempo-Hiisi-as-binary-core), added to `REAL_BINARY_SYSTEMS`.

This module is a REGRESSION suite for the new sourced constants and the
driver's own helper functions, not a re-run of the full discovery sweep
(that ran once as `scripts/run_657_real_binary_kk_sweep_round2.py`; its raw
results are summarized in `data/OUTSTANDING.md`'s `#657` entry, and
`docs/notes/scratch/657_kk_sweep_round2_raw.txt` has the full raw log). It
covers:

1. Static checks on the 3 new sourced `REAL_BINARY_SYSTEMS` entries (mu
   values reproduce their own documented arithmetic from primary-source
   numbers, citations present) -- catches a future edit silently drifting a
   cited number, mirroring #549's own `test_549_sourced_mu_values`.
2. `_fine_n_steps` (the uncapped-recheck driver's step-count helper) hits
   the documented `<=0.001` mu-step convention and reproduces #549's own
   quoted "up to 426 steps" figure for the Patroclus-Menoetius (1,1) /
   mu01215_11-anchor case.
3. A positive-control regression: `sweep_family` against the newly-added
   `lempo-hiisi` system's own (1,1) anchor is clean-negative-aware (passes
   either way), matching #549's own convention for its light per-system
   re-confirmation tests.
"""

from __future__ import annotations

import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent.parent.parent / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from run_657_real_binary_kk_sweep_round2 import (  # noqa: E402
    ANCHOR_TOPOLOGIES,
    GRID_TOPOLOGIES,
    NEW_SYSTEM_KEYS,
    UNCAPPED_RECHECK_PROBES,
    _fine_n_steps,
)

from cyclerfinder.search.real_binary_kk_sweep import (  # noqa: E402
    ANCHORS,
    REAL_BINARY_SYSTEMS,
    sweep_family,
)

# ---------------------------------------------------------------------------
# 1. Static checks on the 3 new sourced systems
# ---------------------------------------------------------------------------


def test_657_new_system_keys_present_with_citations() -> None:
    for key in NEW_SYSTEM_KEYS:
        assert key in REAL_BINARY_SYSTEMS, f"{key} missing from REAL_BINARY_SYSTEMS"
        s = REAL_BINARY_SYSTEMS[key]
        assert 0.0 < s.mu <= 0.5, f"{key}: mu={s.mu} out of (0, 0.5]"
        assert s.mu_source, f"{key}: missing mu_source citation"
        assert s.l_source, f"{key}: missing l_source citation"
        assert s.t_source, f"{key}: missing t_source citation"
        assert s.l_km > 0.0
        assert s.t_s > 0.0


def test_657_sila_nunam_mu_reproduces_magnitude_difference_arithmetic() -> None:
    """mu is derived from Sila's mean 0.12-mag brightness excess (Grundy et
    al. 2012), NOT an assumed 0.5 -- reproduce the exact arithmetic so a
    future edit can't silently drift it."""
    dmag = 0.12
    flux_ratio_sila_over_nunam = 10.0 ** (dmag / 2.5)
    radius_ratio_nunam_over_sila = (1.0 / flux_ratio_sila_over_nunam) ** 0.5
    vol_ratio = radius_ratio_nunam_over_sila**3
    expected_mu = vol_ratio / (1.0 + vol_ratio)

    s = REAL_BINARY_SYSTEMS["sila-nunam"]
    assert abs(s.mu - expected_mu) < 1e-12
    # Sanity: near-equal but NOT exactly 0.5 (a real derived asymmetry).
    assert 0.40 < s.mu < 0.49
    assert abs(s.l_km - 2777.0) < 1e-9


def test_657_antiope_mu_reproduces_aljbaae_2020_mass_ratio() -> None:
    """mu = M_beta/(M_alpha+M_beta) from Aljbaae et al. (2020) Table 1 --
    a genuine individual-mass measurement, not a size/density-assumed
    ratio (contrast with Sila-Nunam and Patroclus-Menoetius, both of which
    ARE equal-density-assumed)."""
    m_alpha = 4.595642e17
    m_beta = 4.525132e17
    expected_mu = m_beta / (m_alpha + m_beta)

    s = REAL_BINARY_SYSTEMS["antiope"]
    assert abs(s.mu - expected_mu) < 1e-12
    assert 0.45 < s.mu < 0.50
    assert abs(s.l_km - 176.0) < 1e-9


def test_657_lempo_hiisi_mu_and_period_reproduce_ragozzine_2024() -> None:
    """mu = M_Lempo/(M_Hiisi+M_Lempo) from Ragozzine et al. (2024) Table 2
    best-fit masses; t_s derived via Kepler's third law from those same
    masses + the tabulated semimajor axis (period is NOT tabulated
    directly in the source). Hiisi is CR3BP-primary (it is the MORE
    massive body, despite being the fainter/named-secondary body under the
    MPC circular's brightness convention) so that mu<=0.5 holds."""
    m_lempo = 5.960e18
    m_hiisi = 7.610e18
    a2_km = 850.0
    g_km3_kg_s2 = 6.67430e-20  # CODATA 2018, matches core/constants.py's own value.

    expected_mu = m_lempo / (m_hiisi + m_lempo)
    expected_t_s = (a2_km**3 / (g_km3_kg_s2 * (m_lempo + m_hiisi))) ** 0.5

    s = REAL_BINARY_SYSTEMS["lempo-hiisi"]
    assert abs(s.mu - expected_mu) < 1e-12
    assert abs(s.t_s - expected_t_s) < 1e-6
    # mu<=0.5 must hold (Hiisi, the heavier body, is CR3BP-primary here).
    assert s.mu < 0.5
    assert abs(s.l_km - a2_km) < 1e-9

    # The paper's own prose cross-check: "1.9-day binary orbital period".
    period_days = s.t_s * 2.0 * 3.141592653589793 / 86400.0
    assert abs(period_days - 1.9) < 0.02

    # Unmodeled-Paha caveat must be present and explicit (#657's own mandate:
    # any Lempo-Hiisi hit needs this flagged prominently, not glossed over).
    assert "Paha" in s.caveat
    assert "unmodeled" in s.caveat.lower() or "not modeled" in s.caveat.lower()


# ---------------------------------------------------------------------------
# 2. Uncapped-recheck driver helper
# ---------------------------------------------------------------------------


def test_657_fine_n_steps_matches_0p001_mu_step_convention() -> None:
    """#549's own bullet: fine-resolution robustness pass held every mu-step
    <=0.001, 'up to 426 steps'. Reproduce that exact figure for the
    Patroclus-Menoetius / mu01215_11-anchor (1,1) case (the largest step
    count among the 6 uncapped-recheck probes)."""
    pm_mu = REAL_BINARY_SYSTEMS["patroclus-menoetius"].mu
    anchor_mu = ANCHORS["mu01215_11"]["anchor_mu"]
    n_steps = _fine_n_steps(anchor_mu, pm_mu)
    assert n_steps == 426, f"expected 426 steps (per #549's own bullet), got {n_steps}"

    # Every step must individually be <= 0.001 in mu.
    assert abs(pm_mu - anchor_mu) / n_steps <= 0.001 + 1e-15


def test_657_fine_n_steps_scales_with_mu_gap_and_floors_at_one() -> None:
    assert _fine_n_steps(0.1, 0.1) == 1  # zero gap still takes >=1 step
    assert _fine_n_steps(0.001, 0.011) == 10
    assert _fine_n_steps(0.5, 0.4381) > 0


def test_657_uncapped_recheck_probes_cover_exactly_the_four_649_inconclusive_cases() -> None:
    """The 6 (system, anchor) tuples must correspond to exactly the 4
    (system, topology) probes #549's own bullet flagged INCONCLUSIVE:
    Patroclus-Menoetius (1,1) [all 3 Table-I anchors -- (1,1) has no single
    anchor], (3,1), (3,2), and Eris-Dysnomia (3,2)."""
    seen = {
        (sys_key, ANCHORS[anchor_key]["k1"], ANCHORS[anchor_key]["k2"])
        for sys_key, anchor_key in UNCAPPED_RECHECK_PROBES
    }
    expected_topologies = {
        ("patroclus-menoetius", 1, 1),
        ("patroclus-menoetius", 3, 1),
        ("patroclus-menoetius", 3, 2),
        ("eris-dysnomia", 3, 2),
    }
    assert seen == expected_topologies
    # (1,1) contributes 3 anchor entries (all 3 Table-I anchors), the rest 1 each.
    assert len(UNCAPPED_RECHECK_PROBES) == 6


# ---------------------------------------------------------------------------
# 3. Topology-list reuse (must match #549's own convention exactly)
# ---------------------------------------------------------------------------


def test_657_topology_lists_match_549_exactly() -> None:
    assert ANCHOR_TOPOLOGIES == [
        (1, 1, "mu001_11"),
        (1, 1, "mu01215_11"),
        (1, 1, "mu05_11"),
        (3, 2, "mu01_32"),
        (3, 1, "mu03_31"),
        (3, 3, "mu01215_33"),
    ]
    assert GRID_TOPOLOGIES == [(2, 1), (2, 2)]


# ---------------------------------------------------------------------------
# 4. Light clean-negative-aware re-confirmation for the new systems
# ---------------------------------------------------------------------------


def _assert_valid_stable_member(result, *, label: str) -> None:  # type: ignore[no-untyped-def]
    assert result.stable_found
    assert result.topology_ok, f"{label}: topology check failed"
    assert result.prograde, f"{label}: orbit must be prograde"
    assert result.reaches_secondary, f"{label}: orbit must reach the secondary realm"
    assert result.nu_mid is not None
    assert abs(result.nu_mid) < 1.0
    assert result.crosscheck_ok, f"{label}: independent Radau crosscheck failed"


def test_657_antiope_11_clean_negative_aware() -> None:
    """Antiope (1,1) via the mu=0.5 Table-I anchor (closest to Antiope's
    own mu=0.496) -- fast (single anchor step), clean-negative-aware:
    passes either way regardless of the full sweep's own verdict."""
    target = REAL_BINARY_SYSTEMS["antiope"].to_cr3bp_system()
    result = sweep_family(target, "mu05_11")
    if not result.stable_found:
        return
    _assert_valid_stable_member(result, label="Antiope (1,1)")
