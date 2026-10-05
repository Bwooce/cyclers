# Lead log — cyclers

## Waiting on the user
- none

## Dispatches (Sydney time)
- 2026-10-05 16:50 AEDT — `corpus-review-fable` (#938): Fable review of the PDF corpus for new routes to novel orbits. Output: `docs/notes/2026-10-05-938-fable-corpus-novel-paths-review.md`. Read-only otherwise.
- 2026-10-05 16:50 AEDT — `ci-keeper-opus` (#939): get CI green (2 failures on origin/main at 0ebb27eb: test_cr3bp_ks_stm Deprit case, test_second_species_continuation LC-vs-KS), then project consistency audit (OUTSTANDING CURRENT STATE, README/data README counts).
  - 3d46e1d5 (local): Deprit lunar fix (elementwise gamma_4 check). 73a: user chose option B (add Cartesian reference, LC 2e-7, LC-vs-KS 5e-7); relayed.

## Handover notes
- none
- 2026-10-05 17:00 AEDT check: ci-keeper local commits 3d46e1d5, 101354f4 (audit; census 392 verified by lead); option B edit in progress in tests/search/test_second_species_continuation.py. Fable: receipt only, no file yet. Push held until B + verification pass.
- 17:05 AEDT: lead approved same-pattern z-block change for 73a (LC vs Cartesian at 1e-9, LC-vs-KS 1e-8), on condition that an mpmath z-block check passes first. Lead commit 66c90b0f (team files).
- 17:10 AEDT: ci-keeper a649314f (73a option B + z-block), 5de65940 (#940 registered). Full verify started 17:02, ETA 17:40-18:20 AEDT; read v_summary.log then. Lead team-files 66c90b0f + uncommitted coordination hazard note.
- 17:15 AEDT: #938 brief extended (user): DOI priority list for corpus gaps; cross-paper technique x case transfer matrix (CSV docs/notes/2026-10-05-938-technique-case-matrix.csv).
- 17:25 AEDT: cyclers_pdf clone was 87 commits behind (71 indexed PDFs missing); lead fast-forwarded to f6f8b4e, 290 -> 431 files. Fable told to source-check [D] claims on top routes. #938 note ETA ~19:30 AEDT.
- 17:21 AEDT: local pytest 50%, 0 failures so far; ETA ~17:40 AEDT.
