# #1025: the literature gate (v3-A1) on the 42 clean #973 members

Status: PRE-REGISTRATION (secs. 1-2), written and committed before the gate runs on any member.
Results in sec. 3. No catalogue writes; nothing in this note is called novel.

Members: the clean two-working-body members of `docs/notes/2026-10-07-973-two-working-body-k-extension.md`.
- gc k = 4: 10 (sec. 10).
- gc k = 5: 6 (sec. 11).
- gc k = 6: 4 (sec. 12).
- ev k = 4: 11 (sec. 13).
- ev k = 5: 11 (sec. 14).

Gate: `literature_check` v3-A1 (`docs/notes/2026-10-08-1045-literature-gate-v3-preregistration.md`
secs. 1-8), offline corpus backend.

## 1. Method

Driver: `scripts/litcheck_1025_973_members.py`. It imports `scripts/litcheck_942_943_scope.py` and uses
its `sig_of`, `from_gauntlet` and `run` unchanged, so each signature is built exactly as for gc-1 and
ev-A:
- primary: Jupiter (gc) or Sun (ev);
- sequence: the time-ordered encounters of one cycle;
- period: k;
- V_inf: per encounter;
- topology: {"repeated-moon"};
- working_bodies: "two" when each body turns >= 0.05 deg somewhere;
- return types: from the cycle key and the leg flight times (FR, FR-n:m, HR, SY, GEN).

The gate runs once per member. In the same run it also re-checks the gate state, with expected values
from the v3-A1 record (`data/942_943_litcheck_scope.json`):
- GanCal#5 and Hollister 1H: "published";
- gc-1, gc-2, ev-C: "inconclusive";
- ev-A, ev-B: "not-found".

The labels were computed first, without the gate (`--labels-only`), so the expectations below can be
specific:
- All 42 are "two".
- No gc member is strictly alternating.
- Every ev member carries at least one FR-n:m return with n:m other than 1:1 (2:1, 3:2 or 2:3), and
  some also a GEN return.

## 2. Expected statuses (fixed before running)

- **gc, 20 members: "inconclusive", naming the Russell & Strange Ganymede-Callisto anchor through F7.**
  - That anchor's body set contains {Ganymede, Callisto}, and it is excluded on working bodies ("one"
    against "two"). This is the gc-1/gc-2 result.
  - The Campagnola 2019 GCGC anchor (two working bodies, `alternating_scope`) is excluded because no gc
    member alternates.
  - **One named risk, gc6-0.** Its time-ordered sequence is (Ganymede, Callisto, Ganymede): two
    distinct Ganymede flybys and one Callisto flyby per cycle, i.e. the cyclic order G G C. The gate's
    #972 F3 rule reads a final encounter equal to the first as a catalogue "closing repeat" and drops
    it. That would make (G, C) look alternating, and the GCGC anchor might then not exclude it.
    - So gc6-0 may come out "published" (GCGC) or with a GCGC-led "inconclusive".
    - If it does, the pre-registered result stands as recorded, and is reported as a
      signature-convention artefact.
    - A diagnostic re-run with the same cycle started at the Callisto encounter, (Callisto, Ganymede,
      Ganymede), is then reported beside it, labelled as a diagnostic. It does not replace the result.
- **ev, 22 members: "not-found".**
  - The H&M anchor (two working bodies, return types {FR, SY, HR}) is excluded on return types. Every
    ev member has an FR-n:m return outside that set, and some also GEN.
  - F7 does not fire: the working labels agree ("two"), so working bodies are not among the exclusion
    reasons. This is the ev-A/ev-B result.
  - No ev member has an HR return.
- **Checks:** as listed in sec. 1, all as expected.

Stop rule: if a check is not as expected, or a member's status differs from the above (other than the
named gc6-0 risk), the run's output is kept, nothing is re-run, and the lead is told.

## 3. Results

(Filled in after the run.)
