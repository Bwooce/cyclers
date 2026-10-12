# #1075 proposed fixes, NOT applied (project paused 2026-10-12 15:10 AEDT)

Saved by the lead from numerics-drift-opus's working tree at the pause. Each patch is a proposal with its
derivation in `docs/notes/2026-10-12-1075-os-update-numerics-drift.md`; apply with `git apply <patch>` at the
resume and run the full gates (ruff, mypy, tests/search, tests/scripts) before committing.

- `free_return_chain_atan2_turn.patch`: compute the turn angle with atan2(|v_in x v_out|, v_in . v_out) instead
  of arccos of the clipped cosine (arccos cannot resolve angles below sqrt(2 * 2**-52) rad = 1.2e-6 deg, the
  proven cause of `test_degenerate_reduces_to_single_free_return_ellipse` failing after macOS 27.0.1).
- `test_saturn_titan_tertiary_newton_tol.patch`: per-hit Newton tolerance for the tertiary known hit, k * sigma
  with k = 3 and the measured section-gap noise sigma = 6.9e-8 (tol 2e-7); the other hits keep 1e-9; identity,
  ghost-margin, Radau and re-approach bounds unchanged.
