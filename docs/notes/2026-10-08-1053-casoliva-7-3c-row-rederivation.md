# #1053: re-derive casoliva-7-3c-em-cycler-2010 as the exact mirror image of casoliva-7-3b (patch prepared; no catalogue edit)

Task `#1053` (earthmoon-opus, 2026-10-08). Lead ruling, option (a): the 7-3c row is to be the EXACT
mirror (x, y, t) -> (x, -y, -t) of the 7-3b row, at the 7-3b row's C, so that the pair carries
Casoliva's printed relation exactly. Background: `#1049` and the correction in
`docs/notes/2026-10-08-1048-1000-coverage-gaps.md` sec. 7.

## 0. Method (pre-registered by the ruling; deterministic, no search)

1. Take the 7-3b row's planar state (x, y, xdot, ydot) and form its mirror (x, -y, -xdot, ydot).
   The mirrored orbit has the same C and period.
2. Propagate it over the 7-3b row's period (DOP853, 1e-12) to check closure.
3. Re-phase it to its y = 0 crossing nearest Casoliva's printed 7-3c initial condition (the
   Table 3 values with the module's 180-degree frame flip, `table3_seed_state`), and store that
   crossing as `state_nd`. It is the same orbit at a different phase, so the row's IC sits at
   Casoliva's own section point.
4. Compute the stability indices with `planar_stability_index` (k_par, k_perp, and Casoliva's
   k_signed, the `#801` convention; 7-3c is not in the `_K_SIGNED_FORCE_PERP` override set, so
   k_signed = max(|k_par|, |k_perp|)).

Script: `scripts/check_1053_7_3c_mirror.py`. Output: `data/1053_casoliva_7-3c_rederive/rederive.json`.

## 1. Results

| quantity | current row (`#797`) | re-derived (`#1053`) | Casoliva printed 7-3c | re-derived vs printed |
|---|---|---|---|---|
| C | 1.0671969118897233 | 1.068655371747616 (= the 7-3b row's; the stored state's own C agrees to 5e-11) | 1.0687623900 | -1.00e-04 relative (`#797`'s 7-3b value; registry mu) |
| T (TU) | 18.850569160664993 | 18.849632202115544 (= the 7-3b row's) | 18.8495559215 | +4.05e-06 relative |
| state_nd (x, y, xdot, ydot) | 0.946182643867633, 4.34e-05, 0.42786762494145186, -1.5132148719027556 | 0.9462675230681173, 0, 0.42795058568089944, -1.5130957243229706 | flip of (-0.9462620909, 0, -0.4279513671, 1.5130807007) | dx 5.4e-06, dxdot -7.8e-07, dydot -1.5e-05 (absolute) |
| full-period closure | 1.47e-09 | 4.4e-09 (re-phased state); 7.9e-09 (unrotated mirror) | - | - |
| k_par / k_perp | 57.043116 / -2.253139 | 57.356153 / -2.269616 (the same as the 7-3b row's: 57.356154 / -2.269616) | k = 57.3519357463 | k_signed +7.35e-05 relative |
| periselene (km) | derive check 13,169.9 | 13,204.2 | 13,210.3 (r_pM) | -4.7e-04 relative |
| perigee (km) | derive check 15,798.8 | 15,710.5 | 15,705.7 (r_pE) | +3.0e-04 |
| apogee (km) | - | 456,992 | 456,901 (r_aE) | +2.0e-04 |
| tof_days_bounds | 81.85821736399389 | 81.85414864034166 | - | - |

- **The re-derived row reproduces Casoliva's printed 7-3c much more closely than the current row.**
  - The IC agrees to about 1e-5. The current row is 8e-5 off in x0 and 1.5e-3 off in C.
  - C and T are reproduced at exactly the 7-3b row's level, which is the point of the exact mirror.
- **Stability (`#1030` conventions).**
  - In-plane: k_par = 57.36, so unstable.
  - Vertical: k_perp = -2.27, outside [-2, 2], so vertically unstable too.
  - The stored stability_index is k_signed = k_par = 57.356153 (Casoliva's Eq. 8 max rule, as for
    7-3b). The row's "UNSTABLE" wording is correct; `#1030`'s wording issue concerns 7-3a and 1-2e
    only.

## 2. Patch set (relative to HEAD at preparation; `git apply --check` passes for both)

- **`data/1053_casoliva_7-3c_rederive/1053_7-3c_rederive.patch`** (`data/catalogue.yaml`).
  - **casoliva-7-3c-em-cycler-2010** changes:
    - jacobi_constant, period_nd, stability_index, state_nd and tof_days_bounds, with their
      comments;
    - the comments on orbit_class (derived periselene), validation_level (per-quantity numbers)
      and casoliva_perigee_km (derived perigee);
    - the period note's T and days;
    - a "CORRECTION (#1053)" paragraph in notes. It says the row was corrected from a member of
      the mirror branch at C 1.067197 (1.5e-3 off in C, within V1 tolerance), cites the `#899`
      test and this note, and says the older `#797` numbers in notes describe the previous state.
  - **casoliva-7-3b-em-cycler-2010:** a "MIRROR PAIR (#1053)" paragraph in notes, citing the `#899`
    test and this note. No value changes.
  - This patch also replaces the withdrawn `#1049` notes patch. That patch stated a wrong
    relation (see the `#1048` note sec. 7).
- **`data/1053_casoliva_7-3c_rederive/1053_validate_evidence.patch`**
  (`src/cyclerfinder/data/validate.py`). It updates the 7-3c V1 entry in `_LEVEL_EVIDENCE` to
  the re-derived per-quantity numbers. It is a text change only; the key and the level are
  unchanged.
  - The source file belongs to the lead's window (not additive core work). I prepared it so the
    evidence string does not go stale.
- **Checks on the patched catalogue:** jsonschema passes; `validate_catalogue` and
  `validate_schema_invariants` return 0 errors (397 rows).
- **Frozen ratchets:** ids, classes, levels and tiers are unchanged, so no census count should
  move. No test pins the 7-3c row's numeric values; I checked with grep for the old C, k and
  state. The tests that use 7-3c (`test_937_turn_gate_fullmodel`,
  `test_earth_moon_class1_resonant_connections`, `test_second_species_continuation`) read
  Casoliva's PRINTED Table 3 values, not the catalogue row.
  - **Expected: none move.** The full set still runs in the window.
- **Docstring note for the lead:** `src/cyclerfinder/search/earth_moon_class1_resonant_connections.py`
  quotes "7-3c (k=57.0431" from the old row. It is descriptive text, not a test. Updating it is
  optional.

## 3. Applied (lead's window, 2026-10-08 17:08-19:05 AEDT, together with `#1042`)

- **Pre-check:** `git diff data/catalogue.yaml` was empty at 931682f7.
- **Patches**, all applied cleanly with `git apply`:
  - `data/1042_ieg_rows/ieg_rows.patch`;
  - this task's `1053_7-3c_rederive.patch`;
  - this task's `1053_validate_evidence.patch`.
- **Validation:** schema OK (397 rows); `validate_catalogue` and `validate_schema_invariants`
  return 0 errors.
- **Full ratchet set** (`tests/data`, `tests/search`, `tests/scripts`,
  `tests/test_catalogue_rediscovery.py`). It ran in foreground chunks under 8 minutes, at machine
  load 5-32. Every module ended EXIT 0, or EXIT 5 for the modules marked entirely slow and
  deselected by the default `-m 'not slow'`.
  - The chunks that hit their own timeout (124) under load were split and re-run, file by file
    where needed.
  - Each slow file then passed alone:
    - `test_504_pluto_charon_kk_sweep` (219 s);
    - `test_earth_moon_class1_resonant_connections` (73c test 298 s, default_k_range 216 s;
      EXIT 0, 1 XPASS of the known cross-platform item);
    - `test_neptune_triton_resonant_families` (255 s);
    - `test_neptune_triton_resonant_connections` (114 s);
    - `test_saturn_titan_resonant_connections`;
    - `test_variational_crnbp_torus` (151 s);
    - `test_variational_qbcp_torus` (91 s).
  - No FAILED or ERROR line in any log.
  - Logs: `<scratch>/earthmoon-opus/win2/`.
- **No census ratchet moved**, as expected.
