"""AST ratchet (#521 phase 2): every scripts/run_*.py must call preflight_search.

The mandatory gate (src/cyclerfinder/data/preflight.py) only prevents a
repeat of the #513/#516/#520 incidents if scripts actually call it. This test
parses every scripts/run_*.py file's AST and asserts a call to
preflight_search appears somewhere in the module -- catching the "forgot to
call the gate" mistake at test time, not 12 hours into an unbudgeted sweep.

_LEGACY_EXEMPT lists every scripts/run_*.py that existed BEFORE this gate was
built (2026-07-02) -- retrofitting all 44 of them is out of scope for #521
phase 2. The four most recent (#515-517, #520 -- the scripts whose incidents
motivated this gate) are NOT exempt: they were retrofitted as the reference
example for how a script should call preflight_search. Any NEW script added
after this test lands must call preflight_search or be added to
_LEGACY_EXEMPT with a comment explaining why -- exemptions are a visible,
reviewed decision, never a silent gap.
"""

from __future__ import annotations

import ast
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_SCRIPTS_DIR = _REPO_ROOT / "scripts"

# Every scripts/run_*.py that predates #521 phase 2 (2026-07-02). Frozen list,
# not a glob-based cutoff, so adding a new legacy-style script always fails
# this test until it's either retrofitted or deliberately added here.
_LEGACY_EXEMPT: frozenset[str] = frozenset(
    {
        "run_299_bifurcation_track_3d_family.py",
        "run_299_lit_check_3d_family.py",
        "run_301_subfamily_validation.py",
        "run_302_aldrin_precursor.py",
        "run_302_s1l1_precursor.py",
        "run_303_bcr4bp_l1_continuation.py",
        "run_303_catalogue_validation_probe.py",
        "run_304_bcr4bp_halo_continuation.py",
        "run_304_catalogue_validation_probe.py",
        "run_305_bcr4bp_gauntlet_probe.py",
        "run_306_phase1_silver_327_gauntlet.py",
        "run_307_precursor_multirev.py",
        "run_320_qp_tori_v3_gauntlet.py",
        "run_325_topology_audit.py",
        "run_330_silver_moontour_v2.py",
        "run_331_silver_v3_gauntlet.py",
        "run_332_silver_v4_gauntlet.py",
        "run_333_qp_family.py",
        "run_335_silver_v41_gauntlet.py",
        "run_338_parallel_demo.py",
        "run_338_silver_v4strict_annual_sweep.py",
        "run_344_stage_b_v2_moontour.py",
        "run_365_mcconaghy_2006_em_k2_v1.py",
        "run_365_russell_ocampo_v1.py",
        "run_392_floquet_low_amplitude.py",
        "run_392_v1_v2_gauntlet.py",
        "run_392_v3_verify.py",
        "run_392_v4_annual_sweep.py",
        "run_392_v4_verify.py",
        "run_393_v1_v2_gauntlet.py",
        "run_393_v3_verify.py",
        "run_393_v4_annual_sweep.py",
        "run_393_v4_verify.py",
        "run_405_cross_system_search.py",
        "run_430_global_precursor.py",
        "run_432_er3bp_discovery.py",
        "run_435_high_e_er3bp.py",
        "run_436_direct_er3bp.py",
        "run_448_region_c_highE.py",
        "run_466_energy_walk.py",
        # #606: a positive-control + pilot DEMONSTRATION of a new capability
        # (variational_periodic_orbit.py), not a catalogue-region discovery
        # sweep -- it has no region_id/n_points to preflight (it re-derives
        # ONE already-known EM L1 halo family member and continues through a
        # documented wall), so the gate does not apply. Visible, reviewed
        # exemption per this file's own docstring.
        "run_606_variational_pilot.py",
        # #608: a bounded generative-seed-model PROOF OF CONCEPT (train on the
        # existing #210 outcome-log corpus, generate N=100 candidates, refine
        # with the existing corrector, compare vs. a uniform-random baseline)
        # -- same category as #606 above: a fixed-N capability demonstration
        # with no region_id/n_points to preflight, not a catalogue-region
        # discovery sweep. Visible, reviewed exemption per this file's own
        # docstring.
        "run_608_generative_seed_poc.py",
        # #614: a bounded family-tagging + nonlinear-encoder COMPARISON against
        # #608's existing baseline (same corpus/split/evaluation pipeline) --
        # same category as #608 above: a fixed-N capability comparison with no
        # region_id/n_points to preflight, not a catalogue-region discovery
        # sweep. Visible, reviewed exemption per this file's own docstring.
        "run_614_family_and_nonlinear_poc.py",
        # #317: a scoping/capability-negative INVESTIGATION (does #317's own
        # real corpus + real sweep-summary records show a learned pre-filter
        # is warranted?), not a catalogue-region discovery sweep -- same
        # category as #606/#608 above: it re-reads existing outcome logs and
        # existing sweep-summary files, with no region_id/n_points to
        # preflight. Visible, reviewed exemption per this file's own
        # docstring.
        "run_317_prefilter_scoping.py",
    }
)


def _calls_preflight_search(tree: ast.Module) -> bool:
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if isinstance(func, ast.Name) and func.id == "preflight_search":
            return True
        if isinstance(func, ast.Attribute) and func.attr == "preflight_search":
            return True
    return False


def _run_scripts() -> list[Path]:
    return sorted(_SCRIPTS_DIR.glob("run_*.py"))


def test_legacy_exempt_list_matches_scripts_that_predate_the_gate() -> None:
    """Sanity-check the frozen list against what's actually on disk.

    Every exempted name must exist; this catches typos and stale entries left
    behind after a script is renamed or removed.
    """
    present = {p.name for p in _run_scripts()}
    missing = _LEGACY_EXEMPT - present
    assert not missing, (
        f"_LEGACY_EXEMPT names files that no longer exist in scripts/: {missing} "
        f"-- update the frozen list."
    )


def test_non_exempt_scripts_call_preflight_search() -> None:
    non_compliant: list[str] = []
    for path in _run_scripts():
        if path.name in _LEGACY_EXEMPT:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if not _calls_preflight_search(tree):
            non_compliant.append(path.name)
    assert not non_compliant, (
        f"the following scripts/run_*.py call no preflight_search() and are not "
        f"in _LEGACY_EXEMPT: {non_compliant}. Either call preflight_search near "
        f"the top of main() (see scripts/run_515_cross_system_3d_search.py for "
        f"the reference pattern) or add the file to _LEGACY_EXEMPT with a "
        f"one-line reason."
    )
