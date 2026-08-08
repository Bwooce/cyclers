# `#797`: Casoliva 2010 Table 3 Earth-Moon resonant-family catalogue writeback

**Task:** `#797`, registered 2026-08-08 during `#787`'s own cross-check against
`data/catalogue.yaml` (found: the Casoliva 2008/2010 / Vaquero 2013 Sec. 4.4 two-body-seeded
p:q-resonant Earth-Moon lineage had ZERO catalogue rows despite `#780` building a full
validation gate module reproducing all 16 Casoliva Table 3 rows). This task writes back the
strongest-sourced members, resolving the two open items the registration bullet explicitly
left unresolved: (a) the admission filter (re-confirmed live, not inherited), and (b) the
`orbit_class` determination (`cycler` vs `quasi_cycler` vs `resonant_po`, checked against each
row's own properties, not assumed).

Source: Casoliva, J., Mondelo, J. M., Villac, B. F., Mease, K. D., Barrabés, E. & Ollé, M.,
"Two Classes of Cycler Trajectories in the Earth-Moon System," *JGCD* 33(5), 2010, pp.
1623-1640, DOI `10.2514/1.46856` ("2010"), Table 3 (p.1630), and its 2008 AIAA 2008-6434
conference precursor. Gate module: `src/cyclerfinder/search/earth_moon_resonant_families.py`
(`#780`), reused unmodified this task.

---

## (a) Admission filter — re-confirmed live from the gate module, not inherited

The registration bullet listed the flagged rows from memory of `#780`'s report. This task
re-ran `table3_row`/`table3_gate_report()` directly (`uv run python -c "..."`, see below) and
confirmed the flags live:

| Excluded (raw filter) | satisfies_resonance | exists_in_em_system | Reason |
|---|---|---|---|
| 1-2a | False | True | fails own p-q resonance relation |
| 1-2b | False | True | fails own p-q resonance relation |
| 2-1c | False | False | fails resonance AND flies through Earth |
| 2-1d | True | False | flies through Earth |
| 3-2a | False | True | fails own p-q resonance relation |
| 3-2d | True | False | flies through Earth |
| 7-3d | False | False | fails resonance AND flies through Earth |

7 of 16 rows excluded by the raw `satisfies_resonance`/`exists_in_em_system` filter (matches
`TABLE3_VALID_DESIGNATIONS`'s 9 survivors: 1-2c, 1-2d, 1-2e, 2-1a, 2-1b, 3-2c, 7-3a, 7-3b,
7-3c). Live confirmation this task ran:

```
$ uv run python -c "from cyclerfinder.search.earth_moon_resonant_families import TABLE3_ROWS; \
    [print(r.designation, r.satisfies_resonance, r.exists_in_em_system) for r in TABLE3_ROWS]"
```

**A stricter bar applied on top of the raw filter, per the task's own standing discipline**
("the row to actually be verified self-consistent (IC/period/Jacobi/stability all check
out)"): `table3_gate_report()`'s own `passed` field requires x0/period/Jacobi AND the
stability index `k` to ALL reproduce Casoliva's own printed Table 3 values within
`TABLE3_IC_GATE_REL_TOL=1e-2`, plus an independent Radau cross-check. Of the 9 raw-admitted
rows, **2 fail this stronger check**: `1-2e` and `7-3a` (both honest misses on the stability
index `k` only — `#780`'s own module docstring documents that IC/period/Jacobi still
reproduce tightly for these two, `1e-5`-`2e-5` relative, but the recovered `k` is wildly
different in both sign and magnitude, unresolved). **These 2 rows are EXCLUDED from
writeback** — a row whose own printed stability character cannot be honestly reproduced is
not verified self-consistent, even though its orbit identity (IC/period/Jacobi) is
independently confirmed. This is registered as a follow-on, `#801` (below).

**Final admitted set: 7 of 16 rows** — `1-2c`, `1-2d`, `2-1a`, `2-1b`, `3-2c`, `7-3b`, `7-3c`.

## (b) `orbit_class` determination — checked per-row against actual properties

The registration bullet flagged this as unresolved and guessed `quasi_cycler` as "the most
likely fit." Direct investigation (prompted by an `advisor()` review before any writeback)
found this guess was wrong: the schema's `quasi_cycler` class requires either a real-ephemeris
`validity_window` (epoch-locked cyclers of opportunity, e.g. the `#569` Uranian moon-pair
rows) or the narrow schema-v5.2 epoch-free CR3BP KAM-corridor carve-out (a torus
*characterizing the neighborhood of an already-catalogued orbit* — not applicable here, since
these ARE the base orbits, not corridors around them). Neither applies to Casoliva's Table 3
orbits.

The real discriminator is `cycler`'s own definition ("transports between encounters") vs.
`resonant_po`'s ("a resonant/libration periodic orbit... with NO demonstrated transport
utility — it never encounters the secondary", schema v4.9, `#453`) — **whether the orbit
actually comes within the lunar sphere of influence**, not whether it satisfies its own
labelled p-q Earth-resonance (all 7 admitted rows do that by construction).

**Method**: for each of the 7 admitted rows, computed the lunar patched-conic SOI
(`a*(m2/m1)^(2/5)`) and Hill radius (`a*(m2/(3*m1))^(1/3)`) at this project's own registry
`mu=0.01215058439469525`: **SOI = 66,182.9 km, Hill radius = 61,524.1 km** (matching the
`66,183 km` figure `em-cycler-21-3d-spatial-2026`'s own row already used for the identical
(2,1)-resonant-orbit question). Cross-checked TWO independent ways per row: (1) Casoliva's own
printed `r_pM` (periselene) column, and (2) an independent DERIVE re-verification this task —
dense propagation (`scipy solve_ivp`, DOP853, `rtol=atol=1e-12`, 4000-point grid) of the
`#780`-corrected orbit over one full period, taking the minimum distance to the Moon. Both
methods agree to `<5e-4` relative for every row:

| Row | Casoliva r_pM (km) | DERIVE min-dist (km) | vs SOI (66,183) / Hill (61,524) | orbit_class |
|---|---|---|---|---|
| 1-2c | 339,553 | 339,569 | 5.1x SOI | **resonant_po** |
| 1-2d | 603,592 | 603,585 | 9.1x SOI | **resonant_po** |
| 2-1a | 90,471 | 90,388 | 1.37x SOI | **resonant_po** |
| 2-1b | 92,590 | 92,590 | 1.40x SOI | **resonant_po** |
| 3-2c | 84,332 | 84,332 | 1.27x SOI | **resonant_po** |
| 7-3b | 13,210 | 13,269 | 0.20x SOI (WELL INSIDE) | **cycler** |
| 7-3c | 13,210 | 13,170 | 0.20x SOI (WELL INSIDE) | **cycler** |

**Result: only 2 of the 7 admitted rows (7-3b, 7-3c) are genuine transport cyclers** — both
sit at `~13,200` km periselene, well inside both the SOI and Hill radius, a real lunar flyby
every period. The other 5 (`1-2c`, `1-2d`, `2-1a`, `2-1b`, `3-2c`) genuinely satisfy their own
p-q Earth-resonance but their closest lunar approach is `1.27`x to `9.1`x the SOI — they never
actually encounter the Moon, exactly `resonant_po`'s own defining carve-out, and exactly
mirroring `em-cycler-21-3d-spatial-2026`'s own established reasoning for its (2,1)-resonant
member (`122,628` km vs the same `66,183` km SOI). Writing all 7 back as `cycler` (the
registration bullet's naive guess) would have asserted a transport claim that Casoliva's own
printed periselene column contradicts for 5 of them.

**Honesty note on `resonant_po`'s stability**: the class's own descriptive prose says "a
STABLE resonant/libration periodic orbit," but the schema-ENFORCED invariant
(`tests/data/test_schema_v47_orbit_class.py`) is only `epoch_locked=false` /
`n_returns=infinite` — stability is not gated. 2 of the 5 `resonant_po` rows written back here
(`1-2d`, `2-1b`) are UNSTABLE per Casoliva's own printed `k` (reproduced to `1.2e-4`/`7.0e-6`
relative). This is stated explicitly on each row rather than silently omitted or forced to fit
the class's own prose.

## `stability_index` convention — a genuine, stated deviation

This catalogue's other CR3BP rows (`braik-ross-*-cycler-2026`, `em-cycler-21-3d-spatial-2026`)
use the Barden half-period-monodromy `nu` convention, which requires a symmetric
perpendicular x-axis crossing. Casoliva's own printed Table 3 IC is a general (frequently
asymmetric) Poincaré `y=0` crossing — her own text: "we have not used this property [the
perpendicular crossings] in this paper" — so Barden's half-period `G`-matrix factorization
does not apply. Each written-back row's `stability_index` instead uses Casoliva's OWN
full-period Eq. 6-8 convention (`k = max(|kappa_par|, |kappa_perp|)`, `#780`'s
`planar_stability_index`), stated explicitly in the field's own comment and a `data_gaps`
entry — **not directly comparable number-for-number to this catalogue's other Barden-nu
`stability_index` values.**

## Rows written back (7)

All in `data/catalogue.yaml`, `source: literature`, `first_published` = Casoliva 2010 JGCD,
`corroborating_sources` = Casoliva 2008 AIAA precursor, `validation_level: V1` (new
`_LEVEL_EVIDENCE` entries added to `src/cyclerfinder/data/validate.py`, same-model CR3BP
reproduction + independent Radau cross-check, mirroring `braik-ross-c11a-cycler-2026`'s own
V1 precedent):

- `casoliva-1-2c-em-resonant-po-2010` — resonant_po, STABLE (k=1.9996)
- `casoliva-1-2d-em-resonant-po-2010` — resonant_po, UNSTABLE (k=4.8579)
- `casoliva-2-1a-em-resonant-po-2010` — resonant_po, STABLE (k=1.9256)
- `casoliva-2-1b-em-resonant-po-2010` — resonant_po, UNSTABLE (k=2.0374)
- `casoliva-3-2c-em-resonant-po-2010` — resonant_po, STABLE (k=1.8752)
- `casoliva-7-3b-em-cycler-2010` — cycler, UNSTABLE (k=57.3519)
- `casoliva-7-3c-em-cycler-2010` — cycler, UNSTABLE (k=57.3519)

Every field traces to Casoliva's own printed Table 3 values (SOURCED) or is explicitly
DERIVE-tagged with a `#797`-cited provenance comment (`jacobi_constant`, `period_nd`,
`stability_index`, `state_nd`, `tof_days_bounds`, `period_days` conversion). `mass_ratio`
uses this project's own DE440-registry `mu` (per `#780`'s own explicit module policy, NOT
Casoliva's own displayed `mu_EM=0.0121529529`, which differs by `1.95e-4` relative — stated
in the field's own comment, not silently substituted). Each row's own `notes:` block
documents the admission reasoning and orbit_class reasoning inline, not just in this note.

## Excluded rows (9 of 16), and why

- `1-2a`, `1-2b`, `3-2a` — fail Casoliva's own `satisfies_resonance` footnote.
- `2-1c`, `7-3d` — fail BOTH `satisfies_resonance` and `exists_in_em_system`.
- `2-1d`, `3-2d` — fail `exists_in_em_system` (fly through Earth) despite satisfying resonance.
- `1-2e`, `7-3a` — pass the raw admission filter but FAIL the stronger self-consistency gate
  (`table3_gate_report().passed=False`): the stability index `k` does not reproduce Casoliva's
  own printed value (IC/period/Jacobi still match tightly, an honest, unresolved `#780` gap).
  Excluded from writeback pending resolution.

## Out of scope (confirmed, not attempted)

Vaquero 2013's own 2:1/3:1 "Earth-Moon Periodic Cyclers" (Sec. 4.4.7) remain explicitly
out of scope — no digit-grade IC exists anywhere for them (`#787`'s digest), tracked
separately as `#799` (direct CR3BP family-continuation reproduction, spawned by `#798`).

## New task registered

**`#801`** — resolve the `1-2e`/`7-3a` stability-index (`k`) reproduction miss `#780`
documented and left unresolved (both rows reproduce IC/period/Jacobi to `1e-5`-`2e-5`
relative but the recovered `k` is wildly wrong in both sign and magnitude — considered and
not confirmed by `#780`: a return-map/sub-period stability convention, or a different
eigenvalue-pair selection specific to the asymmetric printed IC point). If resolved, these 2
rows become writeback-eligible under the same `#797` admission bar (both would be
`resonant_po`: `1-2e`'s own printed `r_pM` is `~1.67` LU / `~642,000` km, `7-3a`'s own printed
`r_pM` is `~1.41` LU / `~541,000` km, both far outside the lunar SOI — neither would be a
`cycler`). Not dispatched.

## Verification

1. `data/catalogue.yaml` round-trips through `yaml.safe_load` + `jsonschema.validate` against
   `data/catalogue.schema.json` — clean, 390 total rows (was 383).
2. `cyclerfinder.data.validate.validate_catalogue()` — 0 errors (the 7 new `V1` declarations
   are covered by the new `_LEVEL_EVIDENCE` entries added this task).
3. `uv run ruff check .` / `uv run ruff format --check .` — clean.
4. `uv run mypy src tests` — clean (839 source files).
5. `uv run pytest tests/data tests/search -q` — full ratchet, see commit log for result.
