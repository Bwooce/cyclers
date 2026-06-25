# #441 (Phase 2) bridge spike — verdict + the ER3BP period_f trap

**Date:** 2026-06-25. Scoping spike for the isolated-ER3BP-family hunt (#440 Phase 2).
Throwaway probes only; no campaign built. Verdict: **NEEDS-WORK — a corrector
convention bug must be fixed first, then the path is tractable.**

## 1. CR3BP→ER3BP e=0 bridge: works exactly
A #440-Phase-1 e=0 resonant member (CR3BP rotating frame, full period T) re-converges
as an exact e=0 ER3BP member via `correct_er3bp_periodic(ER3BPSystem(mu, e=0),
period_f=T/2, is_half_period_residual=True)` — 0 iterations, residual 1e-14, same
orbit. At e=0 the ER3BP IS the CR3BP, so the bridge is exact.

## 2. THE TRAP (the key durable finding) — period_f=T/2 is WRONG for e>0
The ER3BP is non-autonomous: the perturbation is `1/(1+e·cos f)`, period 2π in true
anomaly f. The half-period mirror symmetry that closes a symmetric orbit holds ONLY
when `period_f` is a **multiple of π** (full period a multiple of 2π = the perturbation
period). A Phase-1 member's `T/2` is generally NOT a multiple of π (3:2's is 1.898·π).

Passing `period_f=T/2` to the corrector/continuator at e>0 produces orbits that the
corrector reports as CONVERGED (corrector residual machine-zero) but which **do NOT
actually close** — the independent Radau closure residual grows monotonically (e=0.0016
→ 2.4e-3 … → **1.44 at e=0.30**). The symmetric residual only zeroes (y, xdot) at the
crossing; it never tests full-orbit closure. **Only the independent Radau check catches
the non-closure.** Continuing on T/2 would manufacture a campaign of false-positive
"families."

**Scope of the bug:** this is a Phase-2-specific trap because Phase-1 seeds carry an
arbitrary CR3BP period T. It does NOT retroactively invalidate #432/#435/#436/#437 —
those used commensurate `period_f ∈ {2π, π}` (multiples of π), which close correctly
(the #437 Broucke fold golden used period_f=π and is valid).

## 3. The correct convention + mandatory gate
- Set `period_f` to the resonant commensurate value: the paper's bifurcation periods
  **T₀ = 4π (3/2), 4π/3 (5/2), 2π (3/1), 2π/3 (4/1), 2π (5/1)** (Antoniadou & Libert
  2018; golden-sourced so families line up with the published bifurcation inventory).
- **Gate EVERY member on `independent_residual < 1e-8`, never `corrector_residual`
  alone.** Add a regression test: T/2-continuation's independent residual diverges
  while commensurate-period continuation stays ~1e-12.

## 4. Connected-family continuation (with the fix): reaches e≥0.40 cleanly
3:1 (125 members), 4:1 (101), 5:1, 5:2 all continue from the e=0 seed to e≈0.40 with
all members closing to machine precision, no fold. These are the CONNECTED families
(they have an e=0 limit) — NOT the isolated ones. 3:2 is awkward (q=2 long-libration;
the circular IC is too far — needs a dedicated/digitized seed).

## 5. Strongest first experiment for an ISOLATED family
Isolated families are definitionally DISCONNECTED — you cannot reach them by continuing
the connected family (smooth to e=0.40, no branch-off). Direct seeding from the circular
IC at high e also fails (geometry too far). The path:
1. **Digitize** the paper's high-(e1,e2) configuration per MMR ((0,π) for 3/2 & 5/2;
   (π,0) for 3/1,4/1,5/1) from the Antoniadou & Libert (e1,e2) figures.
2. Converge with `correct_er3bp_periodic` at `period_f=T₀`, **gated on independent
   closure**.
3. **Reverse-continue toward e=0** with the #437 arclength continuator and confirm the
   family DIES before e=0 (no circular limit) — that death is the isolated-family
   signature. (Hard-scrutinise per the #436 lesson: re-verify the death is real, not a
   continuation artifact.)
The figure digitization is the precision-limiting / MEDIUM-confidence step (as the plan
notes). 

## Disposition
- #441 stays IN-PROGRESS. Next build steps: (a) the convention fix + independent-closure
  gate + regression test (bounded safety hardening), then (b) digitize high-e seeds +
  the reverse-continue isolated-family detection.
- The period_f trap is the durable artifact — recorded here so no future ER3BP
  continuation pass repeats it.

## Update 2026-06-25 — recommendation (3.2) independent-closure gate LANDED
Step (a)'s *gate half* is implemented:
`continue_er3bp_family_in_e_arclength(..., independent_gate: float | None = None)`
now stops the walk at the first member whose `independent_residual` exceeds the
gate (and raises `ContinuationError` if the seed itself fails the gate — the
signature of a bad/non-commensurate `period_f`). `None` preserves legacy
behaviour. Regression test `tests/genome/test_er3bp_independent_gate.py`
(slow) pins the mechanism with numbers measured on the Broucke EM Family-7P
family: a 1e-8 production gate preserves the whole valid family (members close
to ~1e-10), while a 2e-10 gate set inside the family's residual range truncates
the walk and every returned member honours the gate.

STILL OPEN for Phase 2: the **full T/2-convention reproduction test** (continue a
#440 Phase-1 resonant seed on `period_f = T/2` and show the independent residual
DIVERGES to ~1.44 by e=0.30 while the commensurate `period_f = T₀` stays ~1e-12)
— needs the #440 seed machinery — plus the high-(e1,e2) digitization and the
reverse-continue isolated-family detection. These belong in the #441 Phase 2
plan, not this spike doc.

## RECONCILIATION 2026-06-26 — #441 CLOSED (method SUPERSEDED by #457; objective re-homed)
This task's planned approach (pulsating-frame `er3bp_continuation` + DIGITIZE the
paper's high-(e1,e2) figure + reverse-continue to confirm death before e=0) is
**RETIRED** in favor of the strictly better #457 method, and its premise was
partly wrong:
- **Premise correction:** §5 assumed isolated families are definitionally
  disconnected and so unreachable by continuation. #457 DISPROVED this for the
  3/1: reading the arXiv:1805.00288 source showed the isolated stable family is
  the HIGH-e STABLE SEGMENT of the UNSTABLE branch II, reachable by continuation
  from the circular bifurcation point (T=2π) — NO digitization, an EXACT seed. The
  "digitize the figure" path (this doc §5, the MEDIUM-confidence step) is therefore
  unnecessary; the precision wall it would have hit is moot.
- **Objective re-homed (nothing dropped):** the isolated-family detection — the
  whole point of #441 — is now done via #457's paper-frame (Antoniadou Eq.1)
  branch-II continuation + JOINT doubly-symmetric corrector (`correct_doubly_
  symmetric_member` in `core/er3bp_paper_frame.py`):
  - **3/1 (π,0): DONE (#457)** — closed to 4.3e-12, STABLE, matches published
    a1=0.480674/e1=0.659951.
  - **4/1, 5/1 (π,0) + 3/2, 5/2 (0,π): TASKED (#460, broadened)** — same method.
- **Durable #441 deliverable RETAINED:** the `period_f=T/2` TRAP finding (§2) +
  the independent-closure gate (`continue_er3bp_family_in_e_arclength(...,
  independent_gate=...)`) + its regression test (`tests/genome/test_er3bp_
  independent_gate.py`) are landed and stay valid — they protect any future ER3BP
  continuation pass from the false-positive-closure trap.
- The CONNECTED-family e>0 continuation (§4) is likewise covered: #442's converger
  tracked the connected 3/1 to e2=0.90.
So #441 is marked COMPLETED as **superseded, objective re-homed to #457/#460** — not
because every MMR's isolated member is reached (3/2, 5/2, 4/1, 5/1 remain, under
#460), but because this WORKSTREAM/approach is closed and its remaining work lives
in #460 under the better method. Marking it done removes a stale ledger entry; the
open isolated-family work is tracked in #460, not lost here.
