# Descriptor-row "off-family" was a ToF artifact — CORRECTION (2026-06-10)

**Supersedes:** the negative in `2026-06-10-dsm-multiarc-closure-results.md` and
the #177 `OFF-FAMILY-AT-ANCHOR-VINF` verdict for the descriptor-bearing
`russell-ch4` rows.

## What happened

I concluded (#180) the 8 descriptor-bearing off-family rows were "0 promotions,
triple-confirmed off-family" — agreement across self-seeding (#177), the local
Takao DSM corrector, and MBH global basin-hopping. **This was a false consensus.**
All three methods inherited ONE upstream bug, so their agreement was not
independent evidence.

## The bug

The Stage-B real-ephemeris Lambert in the closer used `shape.tof_g_days` — the
**coplanar G-arc branch transit** (derived from the idealized circular crossing at
r = Mars SMA = 1.524 AU) — as its time-of-flight, instead of the row's
**tabulated signature transit** (`invariants.transit_times_days`). Real DE440 Mars
at the rendezvous epoch sits at r ≈ 1.40 AU, so the coplanar ToF is the wrong
flight time for the real intercept and forces a high-energy arc → the emerged Mars
V∞ inflates ~1.6–2.1×. The DSM seed inherited the same branch ToF; MBH minimised
ΔV (not anchor-match) from that seed and so never targeted the anchor.

## The evidence it is an artifact

1. **`g_arc_shape` already reproduces BOTH anchors** (coplanar model), independently
   re-verified this session:

   | row | sourced E / M | g_arc_shape emerged E / M |
   |---|---|---|
   | russell-ch4-9.353Gg2 | 9.35 / 10.52 | 9.27 / **10.67** |
   | russell-ch4-6.44Gg3 | 6.44 / 3.74 | 6.36 / **3.84** |
   | russell-ch4-3.78Gg3 | 3.78 / 4.63 | 3.34 / **4.52** |

   The descriptor → (a, e) → V∞ map is sound; the inflation appears only AFTER the
   Stage-B refinement.

2. **A free (epoch, ToF) Lambert to real DE440 Mars, at ~the signature ToF,
   reproduces BOTH anchors** on all examined rows to 0.1–0.3 km/s (per the #1
   re-dig): 9.353Gg2 9.47/10.51 @ 90 d; 9.94Gg3 10.00/10.57 @ 102 d; 3.78Gg3
   3.71/4.72 @ 181 d; 6.44Gg3 6.70/3.75 @ 267 d; 3.64gGg3 3.71/4.52 @ 185 d. The
   2× evaporates.

## The fix

In the Stage-B closure, use the row's **signature transit**
(`invariants.transit_times_days[0]`) as the Lambert ToF — or open ToF as a free
variable in a joint (epoch, ToF) search bracketed near the signature, selecting
the branch whose departure AND arrival V∞ both fall in the anchor band.

## Status / next

These ~6 `russell-ch4` rows are **reachable, on-family, V3-CANDIDATES** — NOT the
recorded negative. The proper next step (a corrected Spec 2): implement the
signature-ToF fix, close each row on real DE440, gate V1 (lamberthub
izzo+gooding + Kepler reprop, emerged V∞ vs sourced anchor) then V3 (n-body
horizon-TCM over laps), and propose V0→V1/V3 writeback for review. The 204 ocampo
rows remain descriptor-gated (no descriptor to seed — genuinely out of reach).

## Lesson

Independent methods that share a common seed or upstream assumption are NOT
independent — their agreement can be a shared-bug artifact, not corroboration. A
confidently-declared NEGATIVE deserves the same adversarial cross-check as a
confidently-declared positive. See memory `orbit-closure-discipline`.
