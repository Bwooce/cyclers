"""M8-UX Phase 2: `cyclerfinder solve` wires optimise_cell_*.

Idealized E-M smoke is fast and deterministic. VEM ballistic convergence is
M-ED-blocked (Task 2.4) — exposed but skipped until M-ED Phases 3-5 land.

Flake note (task #119, 2026-06-06): this module was reported failing once under
xdist parallelism during the M-3D run. Investigation found NO isolation hazard
here — every test uses capsys for stdout and the handlers under test write only
to tmp_path; cli.main() rebuilds its parser per call and holds no module-level
mutable state. The one-off failure was reproduced and traced to a pytest-xdist
*collection mismatch* ("Different tests were collected between gw3 and gwN"),
not to this file: a concurrently-edited untracked test module under tests/data/
(Forge WIP) was written mid-collection, so workers disagreed on the test set and
xdist aborted the entire run (exit 1), surfacing against whichever test the
reporter named. Resolution is to let in-flight WIP test files settle before a
parallel run; nothing to fix in tests/cli. Do not re-investigate this as a CLI
isolation bug."""

from __future__ import annotations

import json

import pytest

from cyclerfinder.cli import main


def test_solve_em_idealized_emits_result(capsys: pytest.CaptureFixture[str]) -> None:
    code = main(
        [
            "solve",
            "--bodies",
            "E,M",
            "--sequence",
            "E-M-E",
            "--k",
            "2",
            "--revs",
            "0,0",
            "--branch",
            "single,single",
            "--fidelity",
            "idealized",
            "--vinf-cap",
            "7.0",
            "--n-starts",
            "2",
            "--no-de",
            "--format",
            "json",
        ]
    )
    assert code == 0
    res = json.loads(capsys.readouterr().out)
    # COMPUTED outputs, clearly labelled as such in the artifact
    assert {"cell_id", "constraints_satisfied", "computed"} <= set(res)
    assert {"closure_residual_kms", "max_vinf_kms"} <= set(res["computed"])


def test_solve_by_cell_id_round_trips(capsys: pytest.CaptureFixture[str]) -> None:
    code = main(
        [
            "solve",
            "--cell-id",
            "EM|E-M-E|k2|r00|bss",
            "--fidelity",
            "idealized",
            "--n-starts",
            "2",
            "--no-de",
            "--format",
            "json",
        ]
    )
    assert code == 0
    assert json.loads(capsys.readouterr().out)["cell_id"] == "EM|E-M-E|k2|r00|bss"


def test_solve_mode_ballistic_is_attempted_on_closed_cell() -> None:
    """2.4 [M-ED]: --mode ballistic plumbs through; arg validation + attempt only.

    Uses a 2-body closed E-M cell so this does not depend on VEM convergence.
    The convergence assertion is the skipped smoke below.
    """
    code = main(
        [
            "solve",
            "--cell-id",
            "EM|E-M-E|k2|r00|bss",
            "--fidelity",
            "ephemeris",
            "--mode",
            "ballistic",
            "--n-starts",
            "1",
            "--format",
            "json",
        ]
    )
    # The call is attempted and returns an exit code (0 ok, or 5 no-candidate);
    # we assert only that the plumbing did not error out the parser.
    assert code in (0, 5)


@pytest.mark.skip(
    reason="OPEN RESEARCH (was M-ED-BLOCKED): the ballistic mode landed "
    "(kwargs verified against cli.py 2026-06-06) but VEM convergence to the "
    "sourced Jones multisets does not occur — dense 16-core scan (task #110) "
    "floors at ~17.9-18.5 km/s vs sourced 2.4-7.0, zero bend-feasible. The "
    "#120 inclined-rung hunt ran and REFUTED 3D-inclination-only: the same grid "
    "on the inclined-circular backend moved the floor only ~0.1-0.4 km/s (still "
    "zero bend-feasible), DE440 control floors at ~18.2. Un-skip when a corrector "
    "variant reaches the family (multi-arc-per-leg / eccentric flybys, NOT "
    "3D-inclination-only); see data/OUTSTANDING.md M-ED open-research entry."
)
@pytest.mark.slow
def test_solve_vem_ballistic_converges() -> None:  # pragma: no cover - skipped
    code = main(
        [
            "solve",
            "--bodies",
            "V,E,M",
            "--sequence",
            "E-M-E-E-V-E",
            "--k",
            "3",
            "--period-basis",
            "E-M",
            "--fidelity",
            "ephemeris",
            "--mode",
            "ballistic",
            "--priority-date",
            "2032-01-01",
            "--vinf-targets",
            "E=5.65,M=3.05",
            "--format",
            "json",
        ]
    )
    assert code == 0


@pytest.mark.slow
def test_solve_em_ephemeris_maintenance_runs(capsys: pytest.CaptureFixture[str]) -> None:
    """Ephemeris-maintenance fidelity over real DE440; asserts shape, not convergence."""
    code = main(
        [
            "solve",
            "--bodies",
            "E,M",
            "--sequence",
            "E-M-E",
            "--k",
            "2",
            "--revs",
            "0,0",
            "--branch",
            "single,single",
            "--fidelity",
            "ephemeris",
            "--mode",
            "maintenance",
            "--priority-date",
            "2032-01-01",
            "--vinf-targets",
            "E=5.65,M=3.05",
            "--n-starts",
            "1",
            "--format",
            "json",
        ]
    )
    assert code in (0, 5)
    res = json.loads(capsys.readouterr().out)
    assert "computed" in res
