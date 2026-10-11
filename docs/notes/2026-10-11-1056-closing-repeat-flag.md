# #1056 Literature gate: a closing-repeat flag on CandidateSignature — PRE-REGISTRATION (2026-10-11), committed before any code

## Problem

The #972 F3 rule (`_has_consecutive_same_body`) drops a final encounter equal to the first before it
reads the sequence cyclically. That is right for the catalogue convention, which repeats the closing
body (e.g. "Callisto-Ganymede-Callisto-Europa-Callisto").

It is wrong for a signature that lists exactly ONE cycle with no closing repeat. gc6-0 shows this
(key `k6|LGanymede>Ganymede/2h|LGanymede>Callisto/2h|LCallisto>Ganymede/2h`):
- Its cycle has two Ganymede encounters and one Callisto encounter. The generator's `sig_of`
  (`scripts/litcheck_942_943_scope.py`) writes the time-ordered list (Ganymede, Callisto, Ganymede).
- F3 drops the last Ganymede as a "closing repeat", leaving (G, C), which reads as alternating.
  The Campagnola 2019 GCGC anchor (alternating, two working bodies) then matches: "published" as run.
- Read correctly the cycle has two consecutive Ganymede encounters (the G -> G 2h return), so it is
  not alternating. The #1025 diagnostic, started at Callisto, reads "inconclusive" naming R-S G-C
  (`data/1025_litcheck_gc6-0_diagnostic.json`).

## The rule chosen: a signature FIELD, not a rotation in `sig_of`

- New `CandidateSignature.closing_repeat: bool | None = None`:
  - None: the legacy F3 inference (drop a final encounter equal to the first). Every existing caller is
    unchanged (catalogue-convention signatures, `signature_from_review_entry`, the campaigns).
  - False: the sequence is exactly one cycle with NO closing repeat. F3 never drops; the cyclic test
    reads the sequence as it is.
  - True: the sequence ends with a closing repeat; drop it (explicit).
- `sig_of` (the generator signatures of the #942/#943/#1025 litchecks, built from one cycle's
  time-ordered flyby table) sets `closing_repeat=False`.
- Why a field and not a rotated start in `sig_of`: a rotation fixes one caller by hiding the convention
  in it. Any other one-cycle caller could still hit the same trap, and a reader of the signature could
  not tell which convention it used. The field states the convention on the signature itself, keeps
  the legacy default, and lets any caller declare it.

## Expected effects (written before running; all must hold, else STOP and report)

1. `scripts/litcheck_942_943_scope.py`: the 26 controls, the negative and the five candidates keep
   their status. `data/942_943_litcheck_scope.json` stays BYTE-IDENTICAL.
   - Only the alternating exclusion reads the flag, and only the Campagnola GCGC anchor declares it.
   - The Jovian controls (GanCal#5, GanCal#1, ...) are "one"-working and are excluded from GCGC by
     working bodies anyway.
   - gc-1 and gc-2 have consecutive Callisto encounters under either reading.
2. `scripts/litcheck_1025_973_members.py` (all 42 #973 members):
   - gc6-0 changes from "published" (GCGC, as run) to "inconclusive", naming the R-S 2009
     Ganymede-Callisto anchor, matching its diagnostic.
   - Every other member keeps its status (the 19 other gc "inconclusive", the 22 ev "not-found").
3. Pinned tests (tests/search/, < 2 s):
   - (Ganymede, Callisto, Ganymede), working "two", closing_repeat=False: not "published"; no citation
     to Campagnola.
   - The same with closing_repeat=None: the legacy F3 behaviour (the #972 F3 test is unchanged).
   - (G, C, G, C, G) with closing_repeat=True: published via Campagnola (the explicit closing repeat).
4. The literature test files stay green; ruff and mypy are clean.

## Process

Implement; tests; run the two scripts; results as a section here; pathspec commits; the hash in a
later step.
