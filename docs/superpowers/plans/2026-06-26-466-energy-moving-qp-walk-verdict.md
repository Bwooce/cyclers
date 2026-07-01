# #466 Energy-Moving QP-GMOS Continuation — Verdict

**Date:** 2026-06-26
**Status:** Built + run. Report-only — NO catalogue writeback, NO novelty claim.
**Lineage:** #290 → #319/#320 (2 SILVER) → #333 (iso-energetic family) → **#466 (energy-moving walk)**
**Design:** `docs/notes/2026-06-26-466-energy-moving-qp-walk-design.md`
**Code:** `src/cyclerfinder/genome/qp_tori_energy_walk.py` + `tests/genome/test_qp_tori_energy_walk.py`
**Campaign:** `scripts/run_466_energy_walk.py` → `data/family_466_energy_walk.jsonl`

## The question #466 had to answer

#333 built a working QP-torus family continuator but it walked the rotation number on
a near-iso-energetic torus (`C_J` span ~6e-7) and never approached the #320 SILVER
Bracket-2 region at `C_J≈3.032`. The **#290↔#320 connection hypothesis** — does the
QP family at `C_J≈3.128` connect down to the SILVER quasi-cycler region at `≈3.032`? —
stayed OPEN. #466 built an energy-MOVING continuation to settle it.

## What was built

**Parent-family-driven energy continuation.** The #296 (1,1) 3D vertical parent family
(`data/family_296_3d_em_11.jsonl`, 265 members) is itself an energy continuation:
`C_J` is strictly monotone in `step_index`, 2.920…3.148. The #333 seed sits at parent
**step 112 (C_J≈3.12785)**; the #320 SILVER Bracket-2 sits at parent **step 8
(C_J≈3.03196)** — 105 members apart on ONE monotone ladder. A genuine Neimark-Sacker
rotation Floquet pair sits on the unit circle at every member from 112 down to 8. The
walk steps the parent orbit down the ladder and re-converges the QP-torus at each
parent energy with `correct_qp_torus`, seeded from that member's monodromy rotation
pair. The parent member's `C_J` IS the torus energy, so the walk MOVES in energy by
construction.

A direct energy-PIN on the torus mean state was prototyped and **rejected** (it fights
the GMOS rows for the same 6 DOF of `c_0` and is compute-infeasible, >75s/step, fails
even at dC=1e-4). Two correctness fixes were required:

1. **Target-direction sign** — the initial `target_step`/`cj_target` early-stop fired
   on the first member (inverted sign). Fixed so the descent stops only once a member
   lies strictly BEYOND the target.
2. **Rotation-pair selection (decisive)** — a 3D vertical-family member carries TWO
   complex unit-circle pairs: a trivial near-+1 central/energy pair (eigen-angle ~0)
   and the genuine Neimark-Sacker rotation pair (angle ~1.66 rad at the seed).
   Selecting "nearest |λ|=1" grabbed the trivial pair arbitrarily and the seeded torus
   **collapsed to the parent periodic orbit (rho→0)** at the very first lower member.
   Fixed by requiring a non-trivial eigen-angle, recovering the rotation pair at all
   105 steps (rho seed 1.5688 at step 112 → 1.3672 at step 8).

## DID THE WALK MOVE IN ENERGY? — YES, decisively

`scripts/run_466_energy_walk.py --stride 4 --target-step 8` →
`data/family_466_energy_walk.jsonl` (+ the step 10/9/8 tail in
`data/family_466_silver_tail.jsonl`). **30 distinct members, C_J spanning
3.127854 (parent step 112) down to 3.031960 (parent step 8) — delta 0.0876**, five
orders of magnitude above the #333 iso-energetic 6.3e-7 floor. Every member is a
genuine irrational torus (all `is_irrational=True`), every member closes the
corrector-independent off-grid topology check at ~2–7e-6 (≪ V1's 1e-4 floor); worst
GMOS Fourier residual 2.5e-4 (one member at step 107, FD-conditioning floor; off-grid
check there is 1.4e-6 — a genuine torus). `terminated=reached_target`. The walk MOVES
in energy by construction and the family is continuous across the whole descent.

Rotation-number trajectory (the slow coordinate): `freq_ratio` rises from 0.250 at the
seed to a **fold at ~0.416 (C_J≈3.072, parent step 39)** then turns back down to 0.217
at the SILVER — a genuine turning point in the rotation number along the energy walk
(the family is NOT monotone in `rho` even though it is monotone in C_J). No
resonance-lock and no mode-truncation breach occurred over the 30 members.

## CONNECTION-HYPOTHESIS VERDICT (#290 ↔ #320) — CONFIRMED

The QP-torus family **continuously connects** the #333 seed (C_J≈3.12785, the #290
Phase-1 Neimark-Sacker torus) down through 30 members to the **#320 SILVER Bracket-2
at parent step 8, C_J=3.03196**. The torus the walk converges at step 8 has frequency
ratio **−0.2169**, reproducing the #320 SILVER Bracket-2's published ratio of **−0.218**
(`docs/notes/2026-06-17-320-first-quasi-cycler-sweep.md`, "Frequency ratio −0.218").
The two SILVER tori from #320 and the #290/#333 seed therefore **lie on ONE
energy-continuation family** — the #290↔#320 connection hypothesis, left OPEN by #333,
is now ESTABLISHED. The path is the #296 (1,1) 3D vertical parent family's energy
ladder; the QP-torus rides it continuously from 3.128 to 3.032.

(Two correctness fixes were required to reach this — see "What was built". The first
descent attempt stopped at a false boundary C_J≈3.084 because the naive
`round(2π/φ)` k-classification rejected valid higher-order-resonance tori; classifying
by the nearest primitive k-th root recovered the full descent. Honest note: the
boundary at 3.084 was NOT a true bifurcation, it was a seeder-classification artifact.)

## NEW QUASI-CYCLER CANDIDATE? — NO new candidate (and correctly so)

The family connects to the **already-identified #320 SILVER Bracket-2** — it does not
surface a NEW quasi-cycler distinct from the #320 SILVERs. The step-8 family member IS
the #320 SILVER Bracket-2 (matching C_J and frequency ratio), which is itself NOT
self-admitted (it has been flagged-pending since #320, awaiting the QP-torus
`literature_check` signature extension + human adjudication). So #466's outcome is
**outcome (b) sharpened into a connection-positive**: the descent reached the SILVER
region AND showed the two regions are one family, but found no torus outside the
already-known #320 set. NOTHING is self-admitted. Were a future finer/longer walk to
surface an off-ladder torus, it would be flagged for human gauntlet review after
`search/literature_check.py`, never written back here.

## Olikara golden status

NOT in the corpus. `grep -ci olikara docs/notes/CORPUS_INDEX.md` → 0; no
Olikara-Scheeres / Olikara 2016 / Henderson-Howell QP-tori family table PDF is present
in the private paper corpus. Registered as a **standing acquisition gap** (per #333
follow-on + `feedback_never_give_up_reproducing_papers`). The walk ships with
self-consistency goldens (energy genuinely moves — the inverse of the #333 containment
assertion; irrationality; corrector-independent off-grid closure). The energy walk was
NOT blocked on the golden, and no golden value was fabricated.

## Honest characterization of the corrector floor

At lower-energy members the N=2 GMOS Fourier-norm residual sits at the ~1.4e-5
FD-Jacobian conditioning floor — just above V1_qp's permissive 1e-5 Fourier gate (the
documented Phase-1 limit). The corrector-INDEPENDENT off-grid topology check stays
~3e-6 (≪ V1's 1e-4 floor), so every walked member is a genuine invariant torus; the
per-member V1_qp Fourier verdict is reported, not used as a hard walk terminator (which
would falsely reject valid tori at the truncation floor).

## Discipline

NO catalogue writeback. NO novelty claim. Any candidate flagged for human gauntlet
review after `search/literature_check.py`, NEVER self-admitted (a QP-torus is never
"novel" until it clears the published record). Frozen census ratchet untouched.
