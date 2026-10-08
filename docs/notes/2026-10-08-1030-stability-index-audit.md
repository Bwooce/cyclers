# #1030: audit of `orbit_elements.cr3bp.stability_index` (which index is stored, and the "stable" wording)

Task `#1030` (earthmoon-opus, 2026-10-08; registered from `#997`, `d7998ba0`). Read-only for
`data/catalogue.yaml`: the wording edits are prepared as patch files and are applied by the lead,
together with the full ratchets, in the morning window.

Script `scripts/run_1030_stability_audit.py`; output `data/1030_stability_audit/audit.json`.
Patches: `data/1030_stability_audit/1030_casoliva_wording.patch` (in scope) and
`data/1030_stability_audit/1030_vaquero_wording_optional.patch` (wider finding; it applies on top
of the first patch).

## 1. Method

- **Producers of the field.** The field is written by hand at writeback from these sources. No
  script writes the catalogue.
  - `search/earth_moon_resonant_families.py::planar_stability_index` / `StabilityIndex.k_signed`
    (the Casoliva rows, `#780`/`#797`/`#801`): Casoliva's Eq. (8), k = max(|k_par|, |k_perp|),
    with the `_K_SIGNED_FORCE_PERP` override {1-2e, 3-2a, 7-3a}. The override stores k_perp even
    when |k_par| is larger.
  - `search/cr3bp_periodic.py::barden_stability` (Barden half index nu = k_par / 2): Vaquero,
    Braik-Ross and the Ross family rows.
  - `ml/seed_generation.py::stability_index` (spectral radius max |lambda|): the `#796` corridor
    rows.
  - Published values: Ross & Roberts-Tsoukkas 2026's "sp".
  - So the field mixes four conventions. `docs/spec.md:615` states only "<= 1 stable; > 1
    unstable", and that is wrong for the Casoliva k (stable iff |k| < 2).
- **Recomputation.** For each of the 50 rows that store a value, I integrated the full-period
  6 x 6 monodromy of the row's own state at the row's own mu (`core.cr3bp.propagate`,
  `stm_mode="fixed_path"`).
  - Planar rows: k_par = tr(M4) - 2 (in-plane) and k_perp = tr(Mz) (vertical); each is stable
    iff it lies in [-2, 2].
  - Spatial rows: the two nontrivial b = lambda + 1/lambda, and the spectral radius.
  - Rows that close worse than 1e-8 (rounded published states: ross-rt-em-cycler-21 and -31) are
    first re-closed at the stored C with `correct_symmetric_fixed_jacobi`.
  - The stored value is identified as whichever candidate it matches to 1e-3 relative.
- **Wording.** "stable"/"STABLE" (not "unstable") is searched in the row's raw YAML, comments
  included.
- **The convention, from the source.** Casoliva et al. 2010, Sec. II.C, Eqs. (6)-(8): k_par =
  lambda_par + 1/lambda_par, k_perp = lambda_perp + 1/lambda_perp, "we define the stability index
  of an orbit by k = max{|k_par|, |k_perp|} ... For 0 <= k < 2 we have linear stability". Her
  text lists "the stable cyclers are 1-2c, 1-2e, 2-1a, 3-2c, and 7-3a".

## 2. Rows whose stored value is k_perp (the brief's list)

| row | stored | = | k_par (in-plane) | k_perp (vertical) | in-plane | vertical | row says | action |
|---|---|---|---|---|---|---|---|---|
| casoliva-1-2c | +1.9996 | k_perp | -1.029 | +1.9996 | stable | stable (marginal) | STABLE | none (k_perp is the max; wording correct) |
| casoliva-1-2e | +1.9998 | k_perp | **-4.191** | +1.9998 | **flip-unstable** | stable (marginal) | STABLE | **reword (patch)** |
| casoliva-2-1a | +1.9256 | k_perp | +1.283 | +1.926 | stable | stable | STABLE | none |
| casoliva-2-1b | +2.0374 | k_perp | +1.513 | +2.037 | stable | **unstable** | unstable (Casoliva) | none (k_perp is the max; the row does not call it stable) |
| casoliva-3-2c | +1.8753 | k_perp | -1.246 | +1.875 | stable | stable | STABLE | none |
| casoliva-7-3a | -1.2985 | k_perp | **-4.965** | -1.2985 | **flip-unstable** | stable | STABLE | **reword (patch)** |

- 1-2c, 2-1a, 2-1b and 3-2c store k_perp because it IS the larger index. That is Casoliva's Eq. (8).
- Only the two `#801` override rows, 1-2e and 7-3a, store the smaller index. On both, the
  in-plane index exceeds 2 in magnitude.
  - Under Casoliva's own Eq. (8), their k would be 4.191 and 4.965, which means unstable.
  - Her printed k reproduces our k_perp to 3e-8 (1-2e) and 4e-4 (7-3a).
  - The `#997` family continuation is an independent computation on an orbit that matches 7-3a's
    crossing state to 1e-5, and it gives b_h = -4.97 there.
- So, evidence first: for these two rows the printed k agrees with the out-of-plane index only.
  Her "stable" verdict is not reproduced for the in-plane block. I do not know whether her table
  carries the out-of-plane index for these rows by design or by slip. The proposed wording states
  both facts and does not change the stored value. That value reproduces the printed number and
  is pinned by the `#801` tests.
- **Correction to my `#997` messages:** I wrote that casoliva-2-1b is "in-plane flip-unstable".
  That is wrong. k_par = +1.513 is in-plane stable; 2-1b is vertically unstable (k_perp = +2.037).

## 3. Wider findings (not in the brief's criterion)

1. **Vaquero 2:1 rows that store Barden nu = k_par / 2 and say STABLE although vertically
   unstable.**
   - vaquero-21-c198: nu = 0.485, k_perp = 2.069 (nu_z = 1.035).
   - vaquero-21-c246: nu = -0.927, k_perp = 2.442 (nu_z = 1.221).
   - The rows print nu_z themselves but still call the orbit "STABLE".
   - Optional patch: "in-plane stable, vertically unstable".
   - vaquero-21-c247 (nu = -1.009), c266, c254 and c313 make no "stable" claim, or claim
     instability correctly.
   - The planar Braik-Ross and Ross (E-M) rows: k_perp lies in [-2, 2], or slightly outside for
     ross-rt-em-cycler-32 (k_perp = -2.014, vertically weakly unstable), whose wording I have not
     reworded.
2. **ross-rt-em-cycler-21-2025: flag RETRACTED (`#1031`).** The first audit re-closed the row at
   the C of its 10-digit rounded state (3.129389531068003), not at the stored C^stable
   (3.129389531088256). The difference, 2e-11, is wider than the row's whole stable window in C
   (1.8e-11), so the audit judged a neighbouring unstable member (nu = -1.35). Fixed in the
   script. On the re-run, the row gives in-plane k_par = +0.048 and k_perp = -1.140: stable, and
   consistent with the stored nu = 0.050. Details and an optional full-precision state patch are
   in `docs/notes/2026-10-08-1031-ross-21-stability-recheck.md`.
3. **The field mixes conventions** (sec. 1).
   - 20 corridor rows (`#796`: 10 Braik-Ross 3-D, 7 Lyapunov 3-D, 3 planar Braik-Ross) store the
     spectral radius (about 1.0). The Casoliva rows store k (stable < 2). The others store nu
     (stable <= 1). The Ross & Roberts-Tsoukkas 2026 abstract-mu rows store the paper's "sp", a
     Barden-type index. em-cycler-21-3d-spatial stores 0.
   - The spec line covers only some of these.
   - Suggested follow-up, not done: a `stability_convention` field, or one convention re-derived
     for every row, with a ratchet. That is a schema change, for the owner or lead.
4. Not computed: none. All 50 rows had a state, period and mass ratio. Apart from the two rows
   re-closed above (1.1e-4 and 1.5e-3), stored-state closure ranges from 2e-13 to 1.4e-5.

## 4. The patches

- `1030_casoliva_wording.patch` (2 rows, 4 lines): casoliva-7-3a and casoliva-1-2e.
  - Changes the notes' "(STABLE, k=...)" and the comment "her own printed STABLE verdict is
    reproduced".
  - New wording: Casoliva's verdict, then this project's vertical and in-plane indices.
  - Unchanged: stability_index values, source_quotes (they quote her convention correctly) and
    names.
  - Checked with `git apply --check` against the current catalogue.yaml.
- `1030_vaquero_wording_optional.patch` (2 rows): vaquero-21-c198 and -c246. Changes the name's
  "(CR3BP, STABLE," and three notes/comment phrases to "in-plane stable, vertically unstable". It
  applies after the first patch.
- I checked that both patches applied in sequence still validate against
  `data/catalogue.schema.json`. The ratchets (names appear in some frozen censuses) are the lead's
  scheduled run.

## 5. Applied (2026-10-08)

Both patches (Casoliva and the optional Vaquero wording) were applied to `data/catalogue.yaml` in the lead's ratchet window of 2026-10-08, in the same commit as the `#970` row and the `#1031` patch.
