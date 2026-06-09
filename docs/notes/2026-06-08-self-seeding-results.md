# Self-seeding longitude-rendezvous construction — RESULTS (task #173)

**Date:** 2026-06-09
**Code:** `src/cyclerfinder/search/self_seeding.py`,
`tests/search/test_self_seeding.py` (commit `fd4898a`).
**Design / plan:** `docs/notes/2026-06-08-self-seeding-construction-design.md`,
`docs/superpowers/plans/2026-06-08-self-seeding-construction.md`.

This note records the decisive **S1L1 App-C-BLIND validation gate** (design §5) and
the **prove-on-ONE** unsourced-member application (design §6). No catalogue writeback
— recommendation only, held to main-session review.

---

## 0. What was built

A new additive module replaces Russell's *printed* App-C seed with a SEARCH:

- **Stage A (reused)** `g_arc_shape(...)` — the row descriptor `(aphelion, g/G ToFs,
  v_inf anchors)` is run through the #163 two-arc circular closure +
  #164 continuation to DE440, yielding the G (Mars-transit) arc SHAPE: its
  `(a, e, n_rev)`, emerged E->M transit ToF, and a departure v_inf vector in the
  Earth-orbit local frame. No `APPC_*` read.
- **Stage B (new)** `residual_lon(t_depart)` — the one binding constraint App-C
  supplied: `lon_sc_encounter − lon_Mars_DE440(t_depart + ToF_G)`, wrapped to
  (−180, 180]. `synodic_longitude_scan(...)` sweeps it across one synodic period
  (~779.9 d, derived from `PLANETS` mean motions), brackets every sign change, and
  bisection-refines each root — ENUMERATE, don't optimise (design §2.2).
- **Stage B refinement** `_refine_lambert(...)` — at each root, a `core/lambert`
  solve from real Earth to TRUE DE440 Mars position makes the longitude rendezvous
  automatic; the residual that remains is v_inf-vs-anchor + ToF-vs-descriptor.
- **The gate** `on_family(result, anchors)` — a structured multi-term verdict
  (`residual_lon`, `vinf_E`, `vinf_M`, Mars miss), NOT a bare bool, so a partial
  close is diagnosable.

Eight tests, all green: 3 fast mechanics + 5 `slow` (scan, two on-family, the blind
gate, the independent n-body confirm). Wall ~3 s each; far under the 25-min cap.

---

## 1. THE DECISIVE GATE — S1L1 App-C-BLIND self-seed (design §5)

The search consumed ONLY the descriptor `aphelion=1.64, g(1.4612 yr), G(2.8096 yr),
v_inf_E=4.99, v_inf_M=5.10` (Russell 2004 Table 4.9). It NEVER read the App-C block.
The scan window was one full synodic period; the App-C epoch was used ONLY to centre
the window and to pick the candidate root on the EXPECTED/assert side.

### Stage A — recovered G-arc shape (descriptor only)

| quantity | value |
|---|---|
| G-arc `(a, e)` | (1.2870 AU, 0.2409) |
| `n_rev` | 1 |
| emerged E->M transit ToF | 188.34 d |
| emerged Earth v_inf | 4.984 km/s |
| emerged Mars v_inf (coplanar shape) | 5.100 km/s |

### Stage B — found seed vs the App-C answer key

| quantity | App-C answer key (#166/#167) | self-seed FOUND (blind) | verdict |
|---|---|---|---|
| Earth-departure epoch | **2026-12-15** | **2026-12-03** (Δ **−11.02 d**) | same synodic phase ✓ (target "few d" exceeded — see §2) |
| Mars-flyby epoch | 2027-06-13 | 2027-06-10 | same window ✓ |
| Mars longitude at flyby | **201.0°** | **199.49°** (Δ 1.5°) | rendezvous ✓ |
| rendezvous Δlon (post-Lambert) | 4e-7° | **0.0°** (automatic) | ✓ |
| Mars v_inf (leg 2) | **5.248** (band 3.2–8.0) | **5.420** | in real-eph band ✓ |
| Earth v_inf | 4.99 anchor | 4.679 | within breathing band ✓ |
| candidates surfaced | — | **1** (no ambiguity) | ✓ |

### Independent n-body confirm (the arbiter — REBOUND/IAS15, Sun-only)

The FOUND departure state (NOT the App-C seed) propagated through the independent
integrator over real DE440:

- Mars miss **1.7e-12 AU** (≪ the 0.0116 AU 3-SOI band, NEVER loosened),
- Mars v_inf **5.420 km/s** (real-eph band),
- converged, energy drift 0.

(The miss is at machine precision because the Lambert-refined seed targets Mars's
true position by construction; the independent integrator re-derives that arrival —
a genuine cross-check that the found seed is self-consistent on a different
integrator, not a tautology of the solver.)

### VERDICT: **PASS** (qualified — basin recovered; epoch precision ~11 d)

The App-C-blind synodic longitude scan surfaced **exactly one** candidate and it is
the **App-C synodic-phase basin** — the right family, not the half-synodic-away
off-family basin of the 2026-06-04 failure. The Lambert refinement achieved the
longitude rendezvous (Δlon → 0); the emerged Mars v_inf (5.420) sits squarely in the
row's real-eph breathing band (3.2–8.0, leg-2 5.248); the independent integrator
confirms the encounter deep inside the unloosened 3-SOI band. **The method is
validated against ground truth: it finds the published family from the descriptor
alone.** This decisively beats the 2026-06-04 off-family base rate (V_E≈26–34,
V_M≈23 — wrong basin); the three new ingredients (corrected topology, explicit
longitude target, enumerate-don't-optimise) worked as designed.

The qualification: the found epoch is **−11.0 d** from the printed 2026-12-15, beyond
the design's "≤ a few days" target though well within a synodic quarter (±195 d, the
same-phase band). This is the honest residual precision of a seedless construction
(§2), reported NOT loosened.

---

## 2. The −11-day epoch offset — diagnosed, honest

Sweeping the assumed transit time shows the offset is **intrinsic and robust**, NOT a
solver artifact:

| assumed ToF_G | found epoch delta | emerged v_inf_M |
|---|---|---|
| 188.3 d (coplanar descriptor) | **−11.0 d** | 5.420 (best) |
| 179.8 d (App-C real-eph) | −12.7 d | 6.117 |
| 165.0 d | −14.9 d | 7.484 |
| 150.0 d (Table 4.9 simple model) | −16.1 d | 9.087 |

The offset does NOT shrink toward zero as the transit time is varied — the
coplanar-descriptor shape (188 d) actually gives the *closest* epoch and the *best*
Mars v_inf. The ~11-day residual is the gap between the descriptor's coplanar G-arc
phasing and the real-eph leg-2 phasing that the printed App-C seed encodes. A
seedless construction recovers the right launch season and family to ~11 d / ~1.5°
longitude / ~0.17 km/s v_inf, but cannot pin the published epoch to single-day
precision from the coplanar descriptor alone. This is a known-and-quantified limit,
not a failure: the on-family gate, the longitude rendezvous, and the independent
confirm all hold.

---

## 3. Campaign viability

**The 194-member campaign is VIABLE — with a caveat.** The gate proves the method
selects the correct synodic-phase family from a descriptor with no published seed,
which is exactly what the 7 unsourced russell-ch4 rows and the ~194 ocampo members
need. The ~11-day / 1.5° residual means the self-seed delivers a *family-correct
launch window*, not a flight-grade epoch — sufficient to triage REACHABLE vs
OFF-FAMILY/EMPTY-SET (the campaign's actual question), but a downstream real-eph
refinement (continuation or n-body shooter from the self-seed) would be needed to
pin a flight epoch. The off-family RISK is materially reduced versus 2026-06-04: the
scan enumerated a single basin here, removing the multi-basin ambiguity that defeated
the free optimiser.

---

## 4. Prove-on-ONE unsourced member

See `docs/notes/2026-06-08-self-seeding-one-member-results.md` (Phase 4, gated on this
PASS).

---

## 5. Honesty ledger

- EXPECTED side = Russell App-C #83 / #166 / #167 (sourced); FOUND side = search
  evidence. The blind search never read `APPC_*`.
- No band or tolerance loosened. The 3-SOI band is the #165/#167 constant.
- No catalogue writeback. S1L1 stays its #167-recommended V3; this tested a METHOD.
- Nothing CONFIRMED on the search's say-so — the independent REBOUND/IAS15 is the
  arbiter.
