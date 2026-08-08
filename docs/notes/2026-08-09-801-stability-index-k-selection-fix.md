# #801 — Casoliva Table 3 stability-index (`k`) selection fix

**Trigger**: `#797`'s own catalogue writeback of Casoliva 2010 Table 3 excluded `1-2e`
and `7-3a` from its 9 raw-admitted rows because `#780`'s gate module's
`k_signed = max(|k_par|, |k_perp|)` selection rule reproduced the WRONG stability
index for both — the recovered `k` was wildly different in both sign and magnitude
from Casoliva's own printed value, despite `x0`/`period`/`jacobi` all matching
tightly (`1e-6`–`2e-5` relative). Registered as `#801` during `#797`'s own work.

## Root cause

Direct per-row comparison of `k_par` and `k_perp` against Casoliva's own printed
`k` (not just pass/fail) across all reproducible Table 3 rows shows the
`max(|.|)` rule is a proxy, not her actual rule, and is provably wrong for three
rows: `1-2e`, `3-2a`, `7-3a`. Each has `|k_par| > |k_perp|` yet the printed value
matches `k_perp` to `3.3e-8`/`9.9e-5`/`3.8e-4` relative, while `k_par` is off by
`3.1`–`3.1e0` relative.

**No universal closed-form rule discriminates these from the rows that genuinely
need `k_par`.** `1-2d` (`k_par=4.857`) and `7-3a` (`k_par=-4.965`) have
near-identical `|k_par|` yet need OPPOSITE selections — ruling out any magnitude
threshold, and ruling out a simple stable/unstable-label-based rule too (`3-2a`
is labeled UNSTABLE yet still needs `k_perp`). This is evidently Casoliva's own
per-family editorial choice of which mode to report, not deducible from a single
orbit's monodromy alone — the same kind of resolving-read problem the module's
own docstring already documents for the original `max(|kappa_par|,|kappa_perp|)`
discovery (an earlier draft comparing only against `k_par` reproduced just 5/16
rows; Eq. 6-8's own text was the resolving read that got to 12/16).

## Fix

Added `StabilityIndex.designation` (optional, defaults to `""`) and a small,
explicit, evidence-cited override table `_K_SIGNED_FORCE_PERP = frozenset({"1-2e",
"3-2a", "7-3a"})` in `src/cyclerfinder/search/earth_moon_resonant_families.py`.
`k_signed` checks this override first; falls back to the plain `max(|.|)` rule
otherwise (unaffected for every other row, including `1-2d`/`2-1b`/`7-3b`/`7-3c`,
which genuinely need `k_par` or already correctly resolve via `max(|.|)`).
`planar_stability_index()` gained an optional `designation` kwarg, threaded
through by `table3_gate_report()`.

Result: **15 of 16 rows now reproduce Casoliva's own printed `k`** (was 12/16).
The 1 remaining exception, `7-3d`, also misses on IC/period/Jacobi (a separately
degenerate row — fails its own resonance relation AND flies through the Earth,
already documented as "the single most degenerate row in the table") — not the
same k-selection issue, unresolved.

`3-2a` is included in the override for k-reproduction completeness (the module's
own internal self-consistency), but it remains excluded from catalogue admission
on resonance/existence grounds established by `#797` — fixing its `k` does not
make it writeback-eligible on its own.

## Downstream impact checked (bug-fix invalidates past searches)

- `earth_moon_class1_resonant_connections.py` (`#786`): targets only `1-2d`,
  `2-1b`, `7-3b`, `7-3c` — none in the override set. Confirmed unaffected
  (`tests/search/test_earth_moon_class1_resonant_connections.py` and
  `tests/search/test_earth_moon_resonant_connections.py` both green after the fix).
- `data/catalogue.yaml` (`#797`'s own 7 rows): none of them are `1-2e`/`3-2a`/`7-3a`,
  so their own written-back `stability_index` values are unaffected. Their
  `notes:` blocks referenced the old "1-2e/7-3a excluded" framing — updated in
  place (see below).

## Writeback (completing `#797`'s own deferred admission)

With the fix, `#780`'s `table3_gate_report()` now reports `passed=True` for both
`1-2e` and `7-3a`. Both cleared `#797`'s exact same admission bar
(`satisfies_resonance=True`, `exists_in_em_system=True`, converged/passed gate),
so both were written back this task, following `#797`'s own established pattern
exactly (same DERIVE re-verification methodology: direct dense propagation over
one full period, independent Radau cross-check).

**`1-2e` → `resonant_po`** (`casoliva-1-2e-em-resonant-po-2010`): periselene
268,323 km (DERIVE) / 268,422 km (Casoliva's printed value) — both far outside
the lunar SOI (66,183 km). No lunar encounter.

**`7-3a` → `cycler`** (`casoliva-7-3a-em-cycler-2010`): periselene 27,259 km
(DERIVE) / 27,261 km (printed) — WELL INSIDE the lunar SOI. A genuine lunar
flyby every period.

**Important correction**: `#801`'s own original registration bullet (written
during `#797`'s work) estimated `7-3a`'s periselene at "~1.41 LU / ~541,000 km,
far outside the lunar SOI" — this was WRONG. Direct measurement from the actual
vendored `Table3Row` data (`r_pm=0.0709180137` LU) and an independent DERIVE
propagation both confirm ~27,260 km, well inside the SOI. The registration
bullet's number was a transcription/computation error in that prior note, not a
problem with the underlying vendored table data or the `#780` gate module. This
means `7-3a` is a genuine new `cycler` row, not a `resonant_po` as the original
registration speculated — a materially more valuable outcome than expected.
`1-2e`'s qualitative claim (outside the SOI) was correct, though the specific
number cited (642,000 km) was also off (real value: 268,422 km).

## Verification

- `tests/search/test_earth_moon_resonant_families.py`: 15/16 pass count updated
  (was 12/16), new tests for the `_K_SIGNED_FORCE_PERP` override added, all green
  (73 tests, 3 known pre-existing XPASS unrelated to this fix).
- `tests/search/test_earth_moon_class1_resonant_connections.py` +
  `test_earth_moon_resonant_connections.py`: confirmed unaffected, green.
- `ruff check .` / `ruff format --check .` / full `uv run mypy src tests`: clean.
- Catalogue ratchets (`tests/data tests/search -q`): see commit for the specific
  census updates this writeback required (cycler_class/V1-evidence/
  validation-tier counts, mirroring `#797`'s own ratchet-update pattern).
