# Self-seeding prove-on-ONE — unsourced member 6.44Gg3 (task #173, Phase 4)

**Date:** 2026-06-09
**Gated on:** the S1L1 App-C-blind gate PASS
(`docs/notes/2026-06-08-self-seeding-results.md`).
**Code/test:** `src/cyclerfinder/search/self_seeding.py`,
`tests/search/test_self_seeding.py::test_self_seed_one_unsourced_member_6_44gg3`.
**Writeback:** NONE. Recommendation only, held to main-session review.

---

## 0. The member and why it was chosen

`russell-ch4-6.44Gg3` — an **UNSOURCED** row: it has NO Russell Appendix-C real-eph
seed block (only 2 of the 9 russell-ch4 parents do, #170). It is the plan's
recommended near-ballistic pick: Russell 2004 Table 4.13 gives aphel 1.54 AU,
descriptor `g(2.087 yr) + G(4.3191 yr)`, v_inf E 6.44 / M 3.74, the longest E->M
transit in the table (262 d, the row's low-v_inf-at-Mars signature), turn ratio
TR=0.95 (near-ballistic — best V3 odds). It is the #164 continuation companion and
already carries `free_return_arcs[]`. The search consumed ONLY this descriptor.

## 1. Stage A — recovered G-arc shape (descriptor only)

| quantity | value | vs anchor |
|---|---|---|
| G-arc `(a, e)` | (1.2626 AU, 0.2849) | — |
| `n_rev` | 2 | — |
| emerged E->M transit ToF | **130.91 d** | vs Russell 262 d (≈ half) |
| emerged Earth v_inf | 6.360 | anchor 6.44 ✓ |
| emerged Mars v_inf (coplanar shape) | 3.842 | anchor 3.74 ✓ |

The coplanar shape matches both v_inf anchors well — but its emerged transit time is
**131 d, half the row's 262-d signature**. That is the seed of the outcome below.

## 2. Stage B — found seed + independent confirm

The scan surfaced **one** candidate over the synodic window:

| quantity | found (blind) |
|---|---|
| Earth-departure epoch | 2026-12-04 |
| rendezvous Δlon (post-Lambert) | 0.0° (automatic) |
| Earth v_inf | 5.994 (anchor 6.44 — in band) |
| **Mars v_inf** | **10.861** |
| Mars miss (Lambert / n-body) | 0 / 1.1e-12 AU (geometrically in-band) |
| on-family verdict | **False** — the v_inf_M term fails |

The Lambert arc reaches Mars at the right longitude (the independent REBOUND/IAS15
confirms machine-precision arrival inside the 3-SOI band), but it arrives at
**v_inf_M ≈ 10.9 km/s**, nearly 3× the row's 3.74 anchor and far above the real-eph
breathing ceiling (~8 km/s). The longitude term passes, the Earth-v_inf term passes,
the **Mars-v_inf term fails decisively**.

## 3. VERDICT: **OFF-FAMILY** (clean, quantified negative — first-class, design §6)

The descriptor → (a,e) map for 6.44Gg3 hosts a longitude rendezvous with real DE440
Mars at this synodic phase — but NOT at the row's anchor v_inf. The mechanism is
explicit: 6.44Gg3 is the LONG-transit / low-Mars-v_inf family (262 d drops v_inf_M to
3.74); the coplanar G-arc shape the continuation returns has a 131-d transit, which
arrives fast and hot (10.9 km/s). The short-transit arc that satisfies the descriptor
ToF/aphelion in circular-coplanar does NOT reproduce the long, low-energy real-eph
arrival the published row requires. This is exactly the design's anticipated risk
(§2.3: "the descriptor → (a,e) map may not host a real-eph longitude solution at the
anchor v_inf") — a clean EMPTY-SET-at-anchor-v_inf negative for this member at this
phase, NOT a bug and NOT a method failure.

This is the honest three-way outcome:
- **CONFIRMED** — no (the Mars v_inf is 3× the anchor).
- **PARTIAL** — no (the failure is family identity, not a TCM-budget overrun).
- **OFF-FAMILY / EMPTY-SET** — **yes.** Recorded as a method-versioned negative.

## 4. Implications

- The S1L1 PASS shows the method WORKS when the descriptor's coplanar shape carries
  the row's real-eph transit geometry (S1L1's 188-d coplanar shape ≈ its 180-d
  real-eph leg). 6.44Gg3 shows it FAILS — cleanly and diagnosably — when the coplanar
  shape's transit time diverges far from the real-eph value (131 vs 262 d).
- The discriminating next step for long-transit rows would be to seed Stage A with a
  multi-rev / longer-ToF G-arc branch (the `n_rev` lever) that reproduces the 262-d
  transit, then re-run the longitude scan. That is a method extension, out of scope
  for this prove-on-one. Recorded for the negative-results registry (method-versioned:
  "self-seed v1, coplanar-shape transit, single G-arc branch").
- **No writeback.** 6.44Gg3 is NOT recommended for promotion; it is a clean negative.

## 5. Honesty ledger

- The descriptor (Russell Table 4.13) is the only input; no fabricated seed.
- No band/tolerance loosened. The OFF-FAMILY verdict was NOT forced to PARTIAL/CONFIRM.
- The independent REBOUND/IAS15 is the arbiter for the geometric arrival; the
  scientific verdict turns on the sourced v_inf anchor, reported not gated.
- One member only — NOT a batch over the 7 / 194 (design §6, out of scope).
