# `#811`: Vaquero 2013 Sec. 4.4.7 Earth-Moon periodic-cycler family catalogue writeback

**Task:** `#811`, registered 2026-08-09 (spawned by `#799`'s successful reproduction). Writes a
selected, justified subset of `#799`'s 129-member reproduction record
(`data/found/799_vaquero_em_cycler_families/results.json`) into `data/catalogue.yaml`, following
`#797`'s Casoliva-writeback pattern field for field. Source: Vaquero Escribano, T. M.,
*Spacecraft Transfer Trajectory Design Exploiting Resonant Orbits in Multi-Body Environments*,
Ph.D. dissertation, Purdue University, 2013, Sec. 4.4.7 ("Earth-Moon Periodic Cyclers"),
pp.169-172, Figs. 4.43-4.44 (no DOI; corpus md5 `fdcbf871322b87cd1dd3448059cb2596`, see
`docs/notes/CORPUS_INDEX.md`). Reproduction evidence: `#799`
(`docs/notes/2026-08-09-799-vaquero-em-cycler-family-reproduction.md`,
`src/cyclerfinder/search/vaquero_em_cyclers.py`), reused unmodified. This closes the
`#787`-identified Casoliva/Vaquero p:q-resonance-lineage catalogue gap on the Vaquero side
(`#797`/`#801` closed the Casoliva side). REPRODUCTION writeback — nothing novel claimed
(`#780` precedent; `literature_check.py` not in play, the reproduction target IS the published
record).

---

## Members selected (6 of 129), and why

The task's own registration set the minimum: the 4 printed-TOF family endpoints + the 2:1
family's stability-transition region. Exactly that minimum is written back — the endpoints are
the only members anchored by a digit-grade Vaquero print (her four prose endpoint TOFs are the
only digit-grade values she prints for these families), and the transition bracket captures the
one family-structure feature her own criterion 4 ("stable or possess small unstable modes")
makes load-bearing. Interior members beyond the bracket add no additional sourced information —
the full 129-member record remains archived in
`data/found/799_vaquero_em_cycler_families/results.json`.

| id | C | role | printed TOF | ours (rel err) | Barden nu / \|lambda\| | periselene (km) | orbit_class |
|---|---|---|---|---|---|---|---|
| `vaquero-21-c198-em-resonant-po-2013` | 1.98 | 2:1 low-C endpoint | 6.39 d | 6.3987 d (0.14%) | 0.4852 / 1.0 (STABLE) | 86,911 | **resonant_po** |
| `vaquero-21-c246-em-cycler-2013` | 2.46 | 2:1 last stable member | — | — | -0.9274 / 1.0 (STABLE) | 49,700 | **cycler** |
| `vaquero-21-c247-em-cycler-2013` | 2.47 | 2:1 first unstable member | — | — | -1.0088 / 1.142 (UNSTABLE) | 48,516 | **cycler** |
| `vaquero-21-c266-em-cycler-2013` | 2.66 | 2:1 high-C endpoint | 4.91 d | 4.9668 d (1.16%) | -2.3270 / 4.428 (UNSTABLE) | 17,675 | **cycler** |
| `vaquero-31-c254-em-cycler-2013` | 2.54 | 3:1 low-C endpoint | 4.90 d | 4.9050 d (0.10%) | 5.6750 / 11.26 (UNSTABLE) | 33,258 | **cycler** |
| `vaquero-31-c313-em-resonant-po-2013` | 3.13 | 3:1 high-C endpoint | 5.04 d | 5.1320 d (1.83%) | 6.6836 / 13.29 (UNSTABLE) | 66,995 | **resonant_po** (SOI-marginal) |

## `orbit_class` determination — per member, not per family

The registration explicitly flagged this as non-trivial, and it is: **the two families do NOT
share one classification, and neither family is internally uniform.** Discriminator: the lunar
patched-conic SOI (66,182.9 km) and Hill radius (61,524.1 km) at the registry
`mu=0.01215058439469525` — the same figures and the same per-row discipline as
`em-cycler-21-3d-spatial-2026` and `#797`'s Casoliva rows.

Checked TWO independent ways per member (`#797`'s discipline — dense propagation, not a table
lookup): (1) `#799`'s own `member_report` dense DOP853 full-period propagation (100k samples +
parabolic extremum refinement), and (2) **`#811`'s own standalone re-derivation** — an inline
CR3BP EOM (no `cyclerfinder` imports at all) integrated with **Radau** (different integrator
family) at `rtol=atol=1e-12`, 20k-sample circular grid + parabolic refinement
(scratchpad script, results below). The two methods agree to **<0.1 km on every quantity for
all six members**; the standalone run also confirms Jacobi drift `<=2.1e-13` and full-period
state closure `<=3.1e-12` per member.

- **2:1 family**: enters the lunar SOI only for `C >= ~2.30` (the family crosses the SOI
  boundary between C=2.29 at 66,768 km and C=2.30 at 65,903 km). So the C=1.98 endpoint
  (86,911 km = 1.313x SOI, 1.413x Hill) is `resonant_po` — it genuinely never encounters the
  Moon — while C=2.46 (0.751x SOI), C=2.47 (0.733x) and C=2.66 (0.267x) are genuine
  transport `cycler`s.
- **3:1 family**: inside the SOI everywhere EXCEPT its C=3.13 top endpoint — the family exits
  the SOI between C=3.12 (65,802 km) and C=3.13 (66,995 km). C=2.54 (0.503x SOI) is a
  `cycler`; **C=3.13 is the closest call in the writeback: 66,995 km = 1.012x SOI, 1.089x
  Hill — outside both, by only ~812 km (1.2%)**. Under the project's SOI discriminator that
  is `resonant_po`; the row states the marginality explicitly (a grid-resolution-conditional
  boundary call, the only 3:1 member outside the SOI in the 60-member record) rather than
  hiding it. Vaquero's own criterion 2 for this family is an L1-LPO connection (this member's
  approach = 1.155x the Moon-L1 distance), satisfied without SOI entry — her "cycler" label is
  her own usage, not this schema's.

Honesty notes carried on rows (per `#797`'s pattern): `vaquero-31-c313` is UNSTABLE despite
`resonant_po`'s descriptive "stable" prose (only `epoch_locked`/`n_returns` is schema-enforced,
same note as `casoliva-1-2d`); the C=2.46/2.47/2.66 members' perigees exceed Vaquero's own GEO
insertion ceiling (`#799`'s criteria-vs-family-extent caveat, restated on each row's
`perigee_km` comment).

## Provenance split (SOURCED vs DERIVE)

- **SOURCED** (Vaquero's own prints, quoted verbatim in each row's `source_quotes`): the two
  families' C-range endpoints and the four endpoint Earth->Moon TOFs (pp.170-171 prose — the
  sentence spans the printed page break), the criterion-4 stability prose (p.170), and the
  free-transfer prose (pp.171-172, quoted on the transition rows — including `#787`'s flagged
  "2.66 < C < 2.54" typesetting-slip reading). `mass_ratio`/`lunit_km`/`tunit_s` are SOURCED
  to the project registry (Vaquero's own mu for Sec. 4.4 is UNSTATED — stated on every row,
  not silently substituted).
- **DERIVE** (`#799`'s continuation, `#811`-cited comments): every IC (`state_nd`), period,
  stability index (Barden half-period nu — the project's usual convention, directly comparable
  to `braik-ross-*` rows and deliberately UNLIKE the `casoliva-*` rows' Casoliva-Eq.-6-8 `k`;
  possible here because `#799`'s members are perpendicular-crossing symmetric orbits),
  periselene/perigee/apogee, TOF values, and day conversions. Interior members' C values are
  additionally flagged (second `data_gaps` entry) as grid-resolution-conditional (`dC=0.01`) —
  Vaquero prints nothing per-member inside the family interior.
- Endpoint rows' `jacobi_constant` is a DERIVED value at a SOURCED target (fixed-Jacobi
  corrector enforcing her printed endpoint C to machine eps; the ~1e-14 recorded offset is
  grid-float accumulation).

Every row carries a `data_gaps` `derive` entry for `state_nd` (no printed IC exists anywhere —
Fig. 4.44's x0 axis is unlabeled, a genuinely non-digitizable graphical source per `#787`).

## Validation level

`V1` on all six rows — same-model CR3BP family reproduction + independent Radau cross-check,
matching the `braik-ross-c11a-cycler-2026` / `casoliva-*` (`#797`) precedent exactly. Six new
`_LEVEL_EVIDENCE` entries added to `src/cyclerfinder/data/validate.py` with per-member numbers
(crossing residuals 7.4e-15..2.0e-12, radau dJ 2.9e-13..2.1e-12, endpoint TOF errors, and the
`#811` standalone-Radau DERIVE re-verification). NOT V2: the four unstable members mechanically
cannot satisfy V2-ballistic bounded-drift; the two stable members simply have not had a
multi-lap run this task (registered as `#823`, below).

## Ratchet-test updates

- `tests/data/test_cycler_class_census.py`: `NON_KEPLERIAN_IDS` +6 (49->55), census
  `non-keplerian` 49->55.
- `tests/data/test_schema_v45_fields.py`: V1-evidence census +6 (V1=37->43), docstring updated.
- `tests/data/test_validation_tier_census.py`: `unvalidated` 113->119 (same
  no-provenance-tags/orthogonal-axis convention as the Casoliva rows, commented in place).
- `tests/search/test_corpus_doi_coverage.py`: NO change needed — the dissertation has no DOI;
  the rows' only cited DOI (Casoliva 2010 corroborating source, `10.2514/1.46856`) is already
  a KNOWN_CORPUS anchor (`#802`).

## Follow-ups registered

- **`#822`**: compute an actual free-transfer (unstable-to-unstable, same-C heteroclinic
  manifold connection) between the 2:1 and 3:1 families in the `C ∈ [2.54, 2.66]` overlap band
  — Vaquero asserts existence (pp.171-172) but prints no transfer; `#799`'s ICs + the existing
  manifold-connection machinery make this a concrete, bounded reproduction target.
  (`#820`/`#821` were claimed by the concurrent `#813` session while this task ran.)
- **`#823`**: V2-ballistic candidacy run for the two linearly stable rows
  (`vaquero-21-c198`, `vaquero-21-c246`): multi-lap bounded-drift long run per the
  `ross-rt-em-cycler-*` V2 precedent.

## Verification

1. `data/catalogue.yaml` round-trips `yaml.safe_load` + `jsonschema.validate` — clean, 398 rows
   (was 392); ids unique.
2. `validate_catalogue()` (via `tests/data/test_validate_catalogue.py`) — clean; the six V1
   declarations are covered by the new `_LEVEL_EVIDENCE` entries.
3. Full `uv run pytest tests/data tests/search -q` ratchet + `ruff check` + `ruff format
   --check` + full `mypy src tests` — see the commit log for results.
