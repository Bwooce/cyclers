# #428 — criteria/physical parameter audit: sourced vs code

**2026-06-23.** Triggered by the discovery that the Russell 200 km flyby floor was
digested (2026-06-07) but the code constant stayed an unsourced 300 km until #426 —
"digest ≠ adoption". This audit classifies every criteria/physical constant used by the
validation + search machinery as:

- **(a) sourced** — traces to a published value with provenance.
- **(b) convention** — an engineering/spec choice, documented as such (legitimate).
- **(c) UNSOURCED DEFAULT shadowing a digested value** — the bug class.

## Findings

| Constant | Value | Class | Provenance / note |
|---|---|---|---|
| `PLANETS[E/M].safe_alt_km` (flyby floor) | 200 km | **(a)** | Russell 2004 p.165 — **was the (c) bug at 300; fixed #426** |
| other planet flyby floors (V/Me/J/S) | various | (a)/(b) | #426 verified defensible; Mercury 1000→200 (#428, BepiColombo demonstrated); see flyby-altitude mining note |
| moon flyby floors (`satellites.py`) | per-moon | (b) labelled | explicitly "convention, not sourced physics"; Callisto 100→200 sourced (#428) |
| band thresholds | <1 / <10 / <300 m/s | **(a)** | band-definitions note + Russell-Ocampo tiers |
| `RUSSELL_BASIS_CYCLES` | 7 | **(a)** | Russell 7-cycle maintenance basis |
| powered_dsm floor / upper | 300 m/s / 3.5 km/s×7 | (a)/(b) | Russell net Δv floor / engineering plausibility upper |
| `V3_BALLISTIC_BUDGET_MPS` | 120 m/s | **(b)** | spec §14 "generic V3 bar" — **SOFT FLAG: internal spec convention, NOT literature-sourced** (defensible: between the <10 and <300 tiers; used only when dv_band is None) |
| `V3_POWERED_MARGIN` | 1.10 | (b) | #175 — 10% margin over documented budget |
| `DRIFT_TOLERANCE_KM` | 50,000 | (b) | plan §4.3 derivation (~0.02°/lap breathing) |
| `REAL_DRIFT_TOLERANCE_KM` | 200,000 | (b) | 4× M6a idealised, documented |
| `V1_TOLERANCE_MPS` | 1e-3 | (b) | spec §14 Lambert-solver agreement bound (numerical) |
| 3.0 km/s/cycle plausibility bar | 3.0 | (b) | engineering bar; validate.py cites sourced ~2.76–2.91 km/s/cycle family Δv |

## Conclusion

**The criteria machinery is in good shape.** Every non-floor constant is either
literature-sourced (a) or a documented engineering/spec convention (b). The **only
genuine class-(c) bug found anywhere was the flyby floor** (E/M 300→200, fixed #426; the
broader per-body floor review is the #428 flyby-altitude mining). No other unsourced
default silently shadowing a digested value was found.

**One soft flag:** `V3_BALLISTIC_BUDGET_MPS = 120` is an *internal spec* convention, not
a literature value — defensible, but it should be labelled "convention, not sourced" so
it is never mistaken for a published bar. (Not a bug; not shadowing a digested value.)

## Process fix (prevents recurrence)

The root cause was a missing **digest→code reconciliation step**: a digest can record a
sourced parameter value while the code constant diverges, with nothing flagging it. Going
forward, any per-paper digest that records a parameter value (altitude floor, budget,
tolerance) must explicitly reconcile against the live code constant — note "code = X,
source = Y" when they differ. Pairs with `feedback_per_paper_digest_todo` and the
update-docs-proactively rule.
