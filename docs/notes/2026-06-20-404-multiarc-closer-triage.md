# #404 — Multi-arc closer triage: which #365 negatives can a closer actually help?

Date: 2026-06-20. Brainstorming #404 ("multi-arc closer substrate") surfaced two facts
that re-scope it before any build.

## 1. The closer is already designed
`docs/superpowers/specs/2026-06-10-dsm-multiarc-closure-lane-design.md` is a full,
**approved** spec for the DSM η-leg multi-arc closure-and-validation lane, with partial
implementation already in tree (`search/dsm_leg.py`, `search/appc_corrected.py`,
`search/bvp_integral.py` — the latter carries the explicit `#388` routing stub for the
deferred Step-6 sun-commensurate-period constraint). #404 is therefore not greenfield.

## 2. A multi-arc genome needs a per-arc descriptor to seed from
The closure lane re-derives a row as two ballistic arcs joined at a flyby — but it can
only do so for rows that publish a **per-arc free-return descriptor** (`free_return_arcs`:
the g/G arc geometry). Rows in the Russell-Ocampo **`n.m.k` summary format** publish only
the cycle-level summary, so there is nothing to seed: they are **descriptor-gated**
(a publication gap), not infrastructure-gated. (The spec's #177 triage already recorded
this; #404 re-confirmed it against the exact #365 rows.)

## Triage of the four #365 negatives

| row | cycler_class | `free_return_arcs` | verdict |
|---|---|---|---|
| russell-ocampo-3.1.2+1 | multi-arc | 0 | descriptor-gated → closer CANNOT seed |
| russell-ocampo-4.3.1-5 | multi-arc | 0 | descriptor-gated → closer CANNOT seed |
| russell-ocampo-4.5.2-2 | multi-arc | 0 | descriptor-gated → closer CANNOT seed |
| **mcconaghy-2006-em-k2** | multi-arc | **2** (generic, 1.46 + 2.81 yr) | **descriptor-bearing → IN SCOPE** |

So the #388/#404 premise ("a multi-arc closer unblocks the 4 negatives") holds for **only
1 of the 4**. The three ocampo rows are blocked by a missing per-arc descriptor (Russell
2004 detailed tables / dissertation), not by missing closer code — their
`negative_results.yaml` resweep conditions were mis-tagged "Implement multi-arc closure
infrastructure" and have been corrected to "acquire the per-arc descriptor."

## Conclusion / next step
The multi-arc closer's reachable #365 payoff is a **single row**, `mcconaghy-2006-em-k2`,
which is already inside the 2026-06-10 spec's descriptor-bearing scope (the spec bounds the
whole lane to ≤ 8 descriptor-bearing rows). The decision (open, for the user): either
**execute that existing approved spec** (writing-plans → build the lane, payoff = the ≤8
descriptor-bearing rows incl. this one), or **defer** the build as low-payoff and treat the
three ocampo rows as the publication-gated negatives they are. No new closer *design* is
warranted — the design exists.
