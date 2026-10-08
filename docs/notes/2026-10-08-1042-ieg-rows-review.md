# #1042 preparation: review of the IEG catalogue rows after #1040, #1041 and #1043

Status: review and patch only; NO catalogue edit (the lead applies the patch in a ratchet window).
Patch: `data/1042_ieg_rows/ieg_rows.patch`, a unified diff against `data/catalogue.yaml` at commit
5da84b2b (the file's last change). `git apply --check` passes. The patched catalogue passes
`validate_catalogue` (0 errors) and `data/catalogue.schema.json`; the current one does too.
The full ratchet suite was not run (lead's window).

## 1. Row `hernandez-2017-jovian-ieg-triple-family` (catalogue lines 9659-9829)

Current:
- orbit_class `cycler`, epoch_locked false, n_returns `infinite`.
- validation_level: field absent (no level recorded).
- model_assumption `circular-coplanar`, trajectory_regime `ballistic`.
- period: pair Io-Europa, k and years null. The note says the natural period is "a small integer
  multiple of Io's period (~1.77 d)" and that the paper was not accessible.
- vinf_kms_at_encounters: all null.
- data_gaps: three CR3BP not-applicable entries.
- notes: family seed; quotes the abstract "repeat indefinitely".

What #1040, #1041 and #1043 change:
- Model note: the paper's ideal model uses the ideal Ganymede period as its synodic period, while its
  radii are built for a rigid 5.2-deg shift. The configuration does not repeat rigidly (per cycle
  Ganymede 0, Europa -20.5, Io -61.5 deg), so its cyclers are quasi-periodic (#1023 sec. 6.3, #1040).
- Reproduction: the Table 4 EGGIE is reproduced in that model at 0.173 km/s, with leg times within
  0.18 d and a gate pass (#1023).
- Repeatability: one cycle ballistic in its own model. Marched cycle by cycle, the gate fails at
  cycle 2 (Io 55.6 deg, ratio 3.72) and no root exists at cycle 3 (#1041). This conflicts with the
  abstract's "repeat indefinitely" and agrees with the paper's own p.10 statement.
- Continuous gravity: undecided. The pinned-end one-cycle chain folds at sigma = 0.505, with the end
  Europa flybys as the growing defect (#1043). This is formulation-limited (#1039).

Proposed (in the patch):
- `period.note`: one added sentence (Table 4 spans 4 synodic periods, printed 28.22 d; reproduced
  cycle 28.02 d = 4 x the ideal Ganymede period).
- `data_gaps`: a new entry `path: n_returns`, `kind: conflict` (abstract "repeat indefinitely" against
  #1041 one cycle), todo_ref #1042. `n_returns` itself is left `infinite` for the family seed pending
  the owner's decision. Changing it would also force an orbit_class change, since cycler requires
  `infinite`.
- `notes`: one added paragraph with the model, reproduction, repeatability and continuous-gravity
  statements above.

Seen, NOT in the patch (for the lead):
- `first_published.note` still says the AAS number "could not be confirmed". The #482 digest
  (`docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md`) confirms
  AAS 17-608.
- `vinf_kms_at_encounters` and `period` could now be filled from Table 4 (sourced: E 9.12,
  G 7.07, I 8.38 km/s; 4 synodic periods). Those are member values on a family-seed row, so the
  member/family split is a schema decision.
- The row carries no validation_level.

## 2. Row `lynam-longuski-2011-ieg-single-period` (catalogue lines 52266-52293)

Current: orbit_class `cycler`, n_returns `infinite`, validation_level V0, trajectory_regime `powered`,
model `circular-coplanar`, period null (the notes state "Laplace period 7.05 d"), a notes string.

What changes: only the reproduction context. #493 reproduced it in the Hernandez ideal model with
T_Laplace = 7.004 d (sourced 7.055 d, -0.72 %). That model has no rigid repeat (#1040). The
rigid-repeat alternative, 7.105 d, is +0.71 % and also inside the 1 % gate. No repeatability or
continuous-gravity test was run for this row. No validation change.

Proposed (in the patch): one sentence appended to `notes`.

Not reviewed (out of scope): `lynam-longuski-2011-gipeipe`, the sibling row. Its notes share the
"#491 ingestion" suffix; the #1040 sentence would apply to it equally if the lead wants it.

## 3. Refresh (2026-10-08, after #1039 and the lead's #1043 record)

The Hernandez row's continuous-gravity sentence in the patch now says the one-cycle EGGIE CLOSES
in continuous gravity in the ideal model (#1043 under #1039 formulation (b): full-mass V_inf
E 9.03, G 7.04-7.05, I 7.88 km/s; altitudes 1,028-7,726 km; identity at 0.5 % mass 0.002 km/s;
IAS15 1e-6 km). It notes the vector-pin fold at 0.505 as an artefact of that formulation, and says
the real ephemeris is untested: the jup365 lane has no validated published control (#968 rung (b),
GanEur#316, not achieved).

The patch was regenerated against `data/catalogue.yaml` at HEAD; the file's last commit is still
5da84b2b. It is 63 lines. `git apply --check` passes against the current working tree, which also
holds another agent's uncommitted edit elsewhere in the file. The patched copy passes
`validate_catalogue` (0 errors) and the JSON schema. Other hunks are unchanged.

Refresh 2 (after #1044): the Hernandez paragraph's real-ephemeris sentence also cites #1044 (GanCal#1
at 2013 did not start: the open-chain formulation fails across its full-revolution return).
Regenerated against HEAD and validated as before. #1044 does not otherwise touch the IEG rows.
