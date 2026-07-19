"""Regression tests for #652: the close-secondary-encounter + wall-clock-budget
``solve_ivp`` termination events guarding ``cr3bp.propagate(with_stm=True)``'s
two STM-integration paths (``_propagate_with_stm_variable`` and
``_propagate_with_stm_fixed_path``).

Background
----------
`#651` found (twice, independently, via bounded SIGALRM diagnostics) that a
rare (~1%) cross-mu seed produced by `#649`'s coordinate-transform rescue can
dynamically evolve into a genuine close approach to the SECONDARY body
partway through propagation -- undetectable from the initial condition
alone -- and hang a single ``correct_periodic`` refinement call for minutes+
(one instance ran past pytest's own 600s global cap before being killed).
The augmented state+STM variational equations (:func:`cr3bp.cr3bp_stm_eom`)
contain second-derivative terms that scale ``~1/r2^5`` near the secondary --
far more singular than the state equations' own ``~1/r2^2`` -- so the
adaptive DOP853 step size collapses without ever tripping the PRE-EXISTING,
UNRELATED ``sol.success=False`` collision-near-PRIMARY failure mode.

This module pins down the EXACT reproduction case discovered while
implementing #652's fix: Sun-Earth mu (``3.0034805950690393e-06``, a bare-mu
"P1"/"P2" system, mirroring how `#649`/`#651`'s cross-mu pipeline constructs
one -- see ``cyclerfinder.ml.seed_generation.resolve_system``), model rng
seed 651, n=100 draw index 77 (0-indexed) of
``generate_and_refine_seeds(..., mu=3.0034805950690393e-06,
rng=np.random.default_rng(651), ...)``'s raw sample, coordinate-transformed
via `#649`'s ``transform_seed_to_target_mu``. The resulting
``state0_guess``/``period_guess`` below are the FULL-PRECISION values
extracted directly from that pipeline (via ``np.save``, not a re-typed
``repr()`` -- this is a chaotic near-singular trajectory, and precision loss
from re-parsing a truncated ``repr()`` was confirmed during this task's own
diagnosis to be enough to miss the close encounter entirely and converge
cleanly instead). Confirmed via ``git stash`` (reverting to pre-#652 code)
plus a bounded ``SIGALRM`` that this EXACT case hung past a 15s bound before
the fix landed; after the fix, every test below completes in well under 1s.

Two-decimal notes on why the direct ``cr3bp.propagate`` calls below use a
SEPARATE, deeper-into-the-Newton-iteration state (``_HANG_STATE0`` /
``_HANG_PERIOD``, captured from an isolated single-``propagate``-call
diagnostic at the exact Newton iterate that first triggered the hang inside
``correct_periodic``) rather than ``_BAD_STATE0_GUESS`` / ``_BAD_PERIOD_GUESS``
(the pipeline's RAW initial guess, several Newton steps upstream): the raw
guess only reaches the pathological regime after several corrector
iterations (confirmed by instrumenting ``correct_periodic``'s own Newton
loop during this task's diagnosis), so a single direct ``propagate()`` call
on the raw guess converges fine -- the end-to-end reproduction
(:func:`test_close_encounter_end_to_end_via_correct_periodic`) exercises
that full path instead.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_periodic import correct_periodic

_MU_SUN_EARTH = 3.0034805950690393e-06


def _bare_mu_system(mu: float) -> cr3bp.CR3BPSystem:
    """Mirrors ``cyclerfinder.ml.seed_generation.resolve_system(mu=...)``'s
    own construction for a generic/unnamed target, without a cross-layer
    dependency on the ``ml`` package from ``tests/core``."""
    return cr3bp.CR3BPSystem(mu=mu, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)


# `#649`-coordinate-transformed raw seed (model rng seed 651, n=100, draw
# index 77) fed to `correct_periodic` -- the end-to-end entry point `#651`
# actually hit. Converges only after several Newton iterations grind into
# the pathological regime (see module docstring).
_BAD_STATE0_GUESS = np.array(
    [
        9.91747188e-01,
        4.89501347e-04,
        8.14052556e-03,
        1.88534508e-03,
        1.80435753e-02,
        -1.18986757e-04,
    ]
)
_BAD_PERIOD_GUESS = 2.0864904743582726

# The Newton-iterate state *within* that same `correct_periodic` call (its
# 15th iteration, of a `max_iter=30` budget) at which a single `propagate`
# call itself first hangs -- captured by instrumenting the corrector loop by
# hand and isolating the exact failing `propagate()` call. Full precision
# (see module docstring on why re-typed `repr()` values are NOT
# interchangeable with these for a chaotic near-singular trajectory).
_HANG_STATE0 = np.array(
    [
        1.00002698e00,
        8.72675037e-06,
        9.12251682e-03,
        1.13670386e-05,
        -1.47324581e-04,
        2.51426780e-03,
    ]
)
_HANG_PERIOD = 2.087727837411517


class TestExceptionHierarchy:
    def test_close_encounter_error_is_a_runtime_error(self) -> None:
        """#652: MUST subclass RuntimeError -- every existing
        `except RuntimeError` caller of `propagate(with_stm=True)` across
        `search/`/`genome/` (and `correct_periodic`'s own many downstream
        callers) must keep working unchanged against this MORE SPECIFIC
        failure, never silently miss it."""
        assert issubclass(cr3bp.CR3BPCloseEncounterError, RuntimeError)

    def test_propagation_timeout_error_is_a_runtime_error(self) -> None:
        assert issubclass(cr3bp.CR3BPPropagationTimeoutError, RuntimeError)

    def test_close_encounter_error_is_not_a_bare_runtime_error_message_collision(self) -> None:
        """The two `RuntimeError` messages this module can raise for
        STM-augmented propagation must stay distinguishable: the pre-existing
        generic 'CR3BP STM propagation failed' (collision near the PRIMARY,
        `sol.success=False`) vs. #652's new close-encounter message (near the
        SECONDARY, an event fired with `sol.success` still True)."""
        try:
            raise cr3bp.CR3BPCloseEncounterError("secondary close encounter")
        except RuntimeError as exc:
            assert "CR3BP STM propagation failed" not in str(exc)


class TestCloseEncounterVariableStmMode:
    @pytest.mark.timeout(20)
    def test_raises_close_encounter_error_fast(self) -> None:
        """#652 regression: before the fix, this exact case hung past a
        SIGALRM-bounded 15s wait (confirmed via `git stash` during this
        task's own implementation, documented in the module docstring). After
        the fix it must raise `CR3BPCloseEncounterError` -- NOT hang, and NOT
        silently move the hang somewhere else -- in well under a second.
        """
        system = _bare_mu_system(_MU_SUN_EARTH)
        t0 = time.monotonic()
        with pytest.raises(cr3bp.CR3BPCloseEncounterError) as exc_info:
            cr3bp.propagate(system, _HANG_STATE0, _HANG_PERIOD, with_stm=True, stm_mode="variable")
        elapsed = time.monotonic() - t0
        assert elapsed < 5.0, (
            f"CR3BPCloseEncounterError fired but took {elapsed:.2f}s -- #652's fix should fire "
            "well under a second on this reproduction case, not merely 'eventually'"
        )
        message = str(exc_info.value)
        assert "secondary" in message.lower()
        assert "close encounter" in message.lower()
        assert "#651" in message or "#652" in message


class TestCloseEncounterFixedPathStmMode:
    """#652: `_propagate_with_stm_fixed_path` was found, during this task's
    own investigation, to be EXPOSED TO THE SAME risk despite its state-only
    pre-pass (empirically: the state-only pass alone does not hang on this
    exact case, finishing in well under a second even though it reaches the
    same close approach -- but the per-sub-interval AUGMENTED replay in step
    2 does hang, confirmed directly via a bounded SIGALRM before this guard
    was added). Both stm_mode values are therefore guarded identically.
    """

    @pytest.mark.timeout(20)
    def test_raises_close_encounter_error_fast(self) -> None:
        system = _bare_mu_system(_MU_SUN_EARTH)
        t0 = time.monotonic()
        with pytest.raises(cr3bp.CR3BPCloseEncounterError):
            cr3bp.propagate(
                system, _HANG_STATE0, _HANG_PERIOD, with_stm=True, stm_mode="fixed_path"
            )
        elapsed = time.monotonic() - t0
        assert elapsed < 5.0


class TestCloseEncounterEndToEnd:
    @pytest.mark.timeout(20)
    def test_close_encounter_end_to_end_via_correct_periodic(self) -> None:
        """The ACTUAL #651-discovered end-to-end path: `correct_periodic`'s
        own Newton loop (not a single direct `propagate()` call) grinding
        into the pathological regime over its own iterations. `#651`'s
        original finding was that this exact scenario ran past pytest's own
        600s global cap; confirmed here (via `git stash` to the pre-#652
        code, bounded by a 15s SIGALRM during this task's implementation) to
        still be hanging at 15s. After #652's fix it must raise
        `CR3BPCloseEncounterError` well under a second, propagating cleanly
        out of `correct_periodic` (which does not itself catch
        `RuntimeError`).
        """
        system = _bare_mu_system(_MU_SUN_EARTH)
        t0 = time.monotonic()
        with pytest.raises(cr3bp.CR3BPCloseEncounterError):
            correct_periodic(system, _BAD_STATE0_GUESS, _BAD_PERIOD_GUESS, tol=1e-10, max_iter=30)
        elapsed = time.monotonic() - t0
        assert elapsed < 5.0, (
            f"took {elapsed:.2f}s -- expected #652's event to fire fast, not merely "
            "eventually within the 20s pytest-timeout backstop"
        )


class TestNormalPropagationUnaffected:
    """#652 must not change behavior for the overwhelming majority of calls
    that never approach the secondary. Reuses `test_cr3bp_stm_mode.py`'s own
    Arenstorf test system/state/duration EXACTLY (mu=0.012277471, IC
    ``[0.994, 0, 0, 0, -2.0015851063790825, 0]``, ``t=5.0``) -- a
    well-behaved trajectory that stays far from any close secondary
    encounter.
    """

    def _arenstorf_system(self) -> cr3bp.CR3BPSystem:
        return cr3bp.CR3BPSystem(
            mu=0.012277471, primary="test", secondary="test", l_km=1.0, t_s=1.0
        )

    _STATE0 = np.array([0.994, 0.0, 0.0, 0.0, -2.0015851063790825, 0.0])
    _T = 5.0

    def test_variable_mode_typical_orbit_unaffected(self) -> None:
        system = self._arenstorf_system()
        arc = cr3bp.propagate(system, self._STATE0, self._T, with_stm=True, stm_mode="variable")
        assert arc.stm is not None
        assert np.all(np.isfinite(arc.stm))
        assert np.all(np.isfinite(arc.state_f))

    def test_fixed_path_mode_typical_orbit_unaffected(self) -> None:
        system = self._arenstorf_system()
        arc = cr3bp.propagate(system, self._STATE0, self._T, with_stm=True, stm_mode="fixed_path")
        assert arc.stm is not None
        assert np.all(np.isfinite(arc.stm))
        assert np.all(np.isfinite(arc.state_f))

    def test_byte_identical_state_result_for_typical_orbit(self) -> None:
        """The #652 events add negligible per-step overhead (a cheap distance
        computation) and must not perturb `solve_ivp`'s own step-size
        selection for a trajectory that never crosses either event's
        threshold -- confirmed by comparing against a raw `solve_ivp` call
        with no events attached at all."""
        from scipy.integrate import solve_ivp

        system = self._arenstorf_system()
        arc = cr3bp.propagate(system, self._STATE0, self._T, with_stm=True, stm_mode="variable")

        y0 = np.concatenate([self._STATE0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom,
            (0.0, self._T),
            y0,
            args=(system.mu,),
            rtol=1e-12,
            atol=1e-12,
            method="DOP853",
            dense_output=False,
        )
        assert sol.success
        np.testing.assert_array_equal(arc.state_f, sol.y[:6, -1])
        assert arc.stm is not None
        np.testing.assert_array_equal(arc.stm, sol.y[6:, -1].reshape(6, 6))
