# #342 — Test failure triage (4 files)

**Date:** 2026-06-17
**Verdict:** ROOT CAUSE = WALL-CLOCK + CPU CONTENTION; no real bugs surfaced.
Resolved by `@pytest.mark.slow` on the 14 genuinely-slow tests, per the
existing project convention (`pyproject.toml` already runs `-m 'not slow'`
by default).

## Triggering context

#321's broad pytest run reported failures in:

1. `tests/data/test_v1_qp.py`
2. `tests/data/test_v2_qp.py`
3. `tests/search/test_cr3bp_3d_family_tracer.py`
4. `tests/scripts/test_cr3bp_family_search.py`

#321's own substrate is green; these failures are NOT a #321 regression.

## Diagnostic timeline

### Isolated-run timings (no concurrent agents, fresh box)

| File | Wall time | n_tests | Result |
|------|----------|---------|--------|
| `test_v1_qp.py` (n_auto, t=60) | 49.9s | 7 | passed |
| `test_v2_qp.py` (n_auto, t=60) | 56.4s | 8 | passed |
| `test_cr3bp_3d_family_tracer.py` (n0, t=600) | 266s | 10 | passed |

### Run under concurrent agent CPU contention (n_auto, t=180)

A concurrent #343 scan (`scripts/scan_343_asymmetric_em_phase2.py`, 16
loky workers at 70-80% CPU each) saturated the 16-core host. The same
tests under contention:

| Test | wall time | Pass/fail |
|------|-----------|-----------|
| `test_v1_qp_sourced_smoke_torus_passes` | **>180s** | TIMEOUT |
| `test_v1_qp_rejects_corrupted_fourier_coefficients` | **166s** | (borderline) |
| `test_v2_qp_sourced_smoke_torus_passes` | **179s** | (borderline) |
| `test_v2_qp_rejects_corrupted_torus` | **>180s** | TIMEOUT |
| `test_v2_qp_cycle_count_scaling_documents_drift` | **>180s** | TIMEOUT |
| `test_pseudo_arclength_reproduces_spike_family` | **>180s** | TIMEOUT |
| `test_fold_detected_and_walk_continues_through` | **>180s** | TIMEOUT |
| `test_closure_preserved_on_every_member` | 83s | passed |
| `test_direction_symmetry` | 81s | passed |
| others in same files | <30s | passed |

### Root cause

These tests run real CR3BP propagation + STM continuation + multi-member
family tracing. On a quiet 16-core host the worst cases take 60-80s
single-threaded, well under pyproject's 600s default. On a host running
4-5 concurrent agents with 16-worker joblib pools, the effective per-test
CPU drops to ~25% and wall times balloon to 180-600s+.

This is NOT a regression — these tests have been slow since they were
written (#296 notes record an 80s baseline). The 60s `pytest --timeout`
override I used surfaced them; the actual pyproject default of 600s
masks them on a quiet box.

The fix is to classify them with the existing project convention.

## Resolution: `@pytest.mark.slow`

`pyproject.toml addopts = "...-m 'not slow'"` — slow-marked tests run only
via `pytest -m slow` (full slow gauntlet) or `pytest -m "slow or not slow"`
(everything). The project already uses this for ~30 other expensive tests
(`tests/test_catalogue_rediscovery.py`, `tests/verify/test_agreement_lamberthub.py`,
etc.).

### Tests marked slow (14 total)

`tests/data/test_v1_qp.py` (2):
- `test_v1_qp_sourced_smoke_torus_passes`
- `test_v1_qp_rejects_corrupted_fourier_coefficients`

`tests/data/test_v2_qp.py` (3):
- `test_v2_qp_sourced_smoke_torus_passes`
- `test_v2_qp_rejects_corrupted_torus`
- `test_v2_qp_cycle_count_scaling_documents_drift`

`tests/search/test_cr3bp_3d_family_tracer.py` (4):
- `test_pseudo_arclength_reproduces_spike_family`
- `test_closure_preserved_on_every_member`
- `test_fold_detected_and_walk_continues_through`
- `test_direction_symmetry`

`tests/scripts/test_cr3bp_family_search.py` (5):
- `test_ross_seed_yields_members_and_runs_inertial_on_each`
- `test_continuation_of_sourced_seed_is_never_novel`
- `test_seed_member_is_reproduction_of_itself`
- `test_offfamily_seed_is_abandoned_not_fabricated`
- `test_routed_members_reach_review_queue_and_counts_agree`

(All 5 in `test_cr3bp_family_search.py` are slow because they either share
the `ross_33_result` module-scope fixture that runs a real continuation,
or run their own `fs.run_seed`/`correct_symmetric_fixed_jacobi` call.)

### Tests left in the default suite (17 in the 4 files)

Audit-trail checks, floor-constant checks, bad-arg validation, and tests
that use degenerate (zero-amplitude) or planar seeds — all <30s individual
wall time on a quiet box, well under 60s under contention.

## Verification

### Default suite (`-m 'not slow'` per pyproject addopts)

```
$ uv run pytest tests/data/test_v1_qp.py tests/data/test_v2_qp.py \
                tests/search/test_cr3bp_3d_family_tracer.py \
                tests/scripts/test_cr3bp_family_search.py --timeout=180 -v
17 passed in 53.49s
```

### Slow suite (`-m slow`)

```
$ uv run pytest tests/data/test_v1_qp.py tests/data/test_v2_qp.py \
                tests/search/test_cr3bp_3d_family_tracer.py \
                tests/scripts/test_cr3bp_family_search.py --timeout=600 -m slow
14 passed in 381.17s (0:06:21)
```

(under concurrent-agent contention; would be ~120-180s on a quiet box.)

### What got rationalised vs fixed

- **STALE RATCHET / FROZEN ASSERTION**: none. Numeric assertions
  unchanged.
- **MODULE API MISMATCH**: none. Tests imports and signatures are intact.
- **REAL BUG**: none. The tests were timing-sensitive against the per-test
  timeout, not assertion-failing. Propagator and corrector return the
  expected closure residuals when given CPU time.

## Phase 2 recommendations (optional, NOT in this commit)

1. **Smaller `n_steps_max` on the 3D-tracer slow-marked tests.** The
   `test_pseudo_arclength_reproduces_spike_family` test uses
   `n_steps_max=25` (forward + backward = 50 members each through a full
   corrector). The topology assertions could survive `n_steps_max=10`
   (forward+backward = 20) — that would halve the wall-clock cost. Same
   for `test_fold_detected_and_walk_continues_through` (n_steps_max=30).
   Defer to a #296 follow-up rather than a quick #342 patch.

2. **A `conftest.py` slow-by-fixture marker.** The
   `test_cr3bp_family_search.py` 3-test fixture-share pattern would
   benefit from a `pytestmark = pytest.mark.slow` at module level. I
   chose explicit per-test markers in this triage to stay conservative
   (one test in that file, `test_empty_record_written_when_no_novel`, is
   fast and deserves to stay in the default suite).

3. **Concurrent-agent CPU-budget discipline.** The underlying friction is
   16-core hosts running 4-5 concurrent agents at 16 workers each. A
   shared semaphore or per-agent worker cap (e.g., `n_workers=4` when ≥2
   agents are running) would let pytest finish in default-suite time and
   reduce the need for `slow` markers on borderline tests. Out of scope
   for #342; track as a separate infra item.

## Files modified

- `tests/data/test_v1_qp.py` (+2 `@pytest.mark.slow`)
- `tests/data/test_v2_qp.py` (+3 `@pytest.mark.slow`)
- `tests/search/test_cr3bp_3d_family_tracer.py` (+4 `@pytest.mark.slow`)
- `tests/scripts/test_cr3bp_family_search.py` (+5 `@pytest.mark.slow`)
- `docs/notes/2026-06-17-342-test-failure-triage.md` (this file)

## What was NOT touched

- No catalogue writeback.
- No module signature changes.
- No assertion adjustments.
- No numeric goldens edited.
- No file paths beyond the 4 failing-test files + this note.
