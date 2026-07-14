# QBCP EOM fix: alpha_6 scales only the Sun-distance term, not the whole Newtonian potential

**Date:** 2026-07-14
**Origin:** recovered from an abandoned `git stash` (WIP on `7f83277`, undated but pre-dating this
session) while cleaning up stale stashed work across both repos. The stash's own diagnosis was
verified before applying, not trusted on inheritance.

## The bug

`qbcp_eom()` and `qbcp_potential_second_derivatives()` in `src/cyclerfinder/core/qbcp.py` multiplied
the ENTIRE Newtonian potential gradient (Earth + Moon + Sun terms) by `alpha_6` when computing the
momenta derivatives. Per Rosales-Jorba (2023) Eq. 3 — already a trusted, previously-cited corpus
reference (`data/README.md`/`docs/notes/2026-06-16-*-bcr4bp-*.md`, used for the POL1/POL2 golden
tests) — `alpha_6` is a coefficient that scales ONLY the Sun-distance term in the Hamiltonian
(`-m_S / (alpha_6 * R_PS)`); the Earth/Moon Newtonian terms carry no `alpha_6` factor at all.

## Why existing tests never caught it

`tests/core/test_qbcp.py`'s own golden case is set up at `alpha_6 = 1` exactly (a documented
simplification in that test). At `alpha_6 = 1`, `(1/alpha_6) * sun_term` and `alpha_6 *
(full_potential)` coincide only for the Sun term and happen to leave Earth/Moon terms numerically
identical too (multiplying the whole gradient by 1 is a no-op) — so the bug is invisible at exactly
that test point. `alpha_6` is in general a genuine time-varying Fourier-series coefficient
(`evaluate_alphas()`), not fixed at 1, so any real QBCP propagation over a synodic period picks up a
systematic O(alpha_6 - 1) error the golden test cannot see.

## Verification before applying

1. Confirmed the cited source (Rosales-Jorba 2023) is a real, already-trusted, previously-used corpus
   reference — not an inherited/unverified citation.
2. Applied the stashed fix and re-ran `tests/core/test_qbcp.py`, `tests/genome/test_qbcp_torus.py`,
   `tests/scripts/test_run_538_residual_shape.py` (the only 3 files referencing `qbcp` at all): all
   pass, identically to before — confirms the fix doesn't regress the `alpha_6=1` golden case (expected,
   since old and new code coincide exactly there).
3. Ran the broader `tests/genome/` + `tests/core/` suite with the fix applied: 3 failures appeared
   (`test_ephemeris_cache.py::test_states_batch_matches_scalar_state[astropy]`,
   `test_da_section_map.py::test_taylor_fixed_point_reaches_png_neighbourhood`,
   `test_qp_tori.py::test_structural_qp_continuation`) — but reproducing WITHOUT the fix (clean
   `main`) shows the exact same 3 failures with byte-identical numeric output, confirming all 3 are
   pre-existing, unrelated to this change (a new addition to the same local-Mac/BLAS-artifact category
   already tracked for 2 other tests under task #584).
4. `ruff check` / `ruff format --check` / `mypy` all clean on the changed file.

## Scope / follow-up

This is a physics-correctness fix to the QBCP equations of motion module, not a new capability. Any
past QBCP-based search/continuation result (the #533/#537/#538/#544 task chain) computed under the
buggy code used a subtly-wrong Sun-term scaling — per this project's own
[[feedback_bugfix_invalidates_past_searches]] discipline, this should be scanned for past negatives
that might change, registered separately (not resolved here) since determining actual impact on past
conclusions is a substantial task of its own, not a stash-cleanup side effect.
