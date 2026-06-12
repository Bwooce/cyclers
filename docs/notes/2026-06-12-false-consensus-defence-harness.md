# False-consensus defence harness (#202)

One-stop index for the executable false-consensus doctrine. The doctrine
(memory `orbit-closure-discipline`, "FALSE-CONSENSUS DEFENCE" addendum;
`docs/notes/2026-06-11-project-review-results.md` §"The false-consensus
doctrine") turned three real incidents into operational rules. This note links
the artifacts that make it checkable.

**Principle:** agreement between N checks is only worth what they do NOT share.

## The three incidents this harness defends against

| Incident | Shared component that defeated "independence" |
|---|---|
| #180 | three "independent" methods inherited one upstream ToF bug |
| #197 | `crosscheck_leg` read its Lambert endpoints from the artifact under test |
| #198 | a 63 s UTC/TDB epoch-conversion offset shared between primary path and its cross-check |

## Fault-injection suite (doctrine item 3 — the core deliverable)

`tests/verify/test_fault_injection.py`. For each SHARED component a validation
gate depends on, the suite poisons it via `monkeypatch` (source is never
mutated) and asserts a SPECIFIC gate FIRES. Each poison shows BOTH halves of
rejection power — a CLEAN test (gate passes unpoisoned, so it is not vacuous)
and a POISONED test (gate fails on the fault, so it has teeth).

| # | Shared component | Poison | Gate asserted to fire | Incident |
|---|---|---|---|---|
| 1 | Epoch convention (`_J2000_TDB_JD`) | shift J2000(TDB) ref +60 s | `test_epoch_anchor.py` Horizons absolute-epoch anchor | #198 |
| 2 | Frame handedness (`_InclinedCircularBackend._rotation`) | flip `R_x(+inc)` → `R_x(-inc)` | `test_ephemeris_inclined.py` DE440 orbit-normal anchor | #199 |
| 3 | Crosscheck endpoint (shared `(r1, r2)`) | displace embedded encounter 1000 km | `test_crosscheck.py::test_crosscheck_catches_poisoned_endpoint` (REFERENCED, not duplicated; it is the template) | #197 |
| 4 | Earth-Moon μ (registry GMs) | perturb `PRIMARIES["Earth"]` +1% | `test_cr3bp.py::test_earth_moon_mu_physical` (sourced μ = 1.2150584270572e-2) | #212a |
| 5 | Signature / transit ToF | shift a leg's claimed arrival epoch +30 d with stale geometry | crosscheck endpoint-independence gate (`endpoint_mismatch_km`) | #180 |

### Finding (an undefended shared component, recorded not papered over)

`crosscheck_leg(..., independent_endpoints=False)` is an explicitly undefended
shared component: with the endpoint re-query disabled, a 1000 km poison sails
through with solver agreement intact. This is recorded — and the mode is kept
OFF by default — in
`tests/verify/test_crosscheck.py::test_crosscheck_escape_hatch_reproduces_shared_endpoint_blindness`.
No in-tree validation path may turn it off.

## Positive-control convention (doctrine item 4)

Before a *negative* result is trusted, the method must re-find a KNOWN solution
*through the identical pipeline configuration* (the rule that would have killed
#180 — a mis-configured method that finds nothing is uninformative, not a real
negative).

Reusable assert helper: `cyclerfinder.verify.positive_control.assert_positive_control`.
A campaign wraps its solver as a zero-arg closure that runs a known-positive
control with the SAME config object the negative sweep then uses, and asserts
re-discovery first. Sharing the config by construction is the point: a passing
control certifies the pipeline, not just the math. This is opt-in plumbing —
it does not retrofit every campaign.

## Consistency-vs-independence + "shared with primary path:" (doctrine items 1-2)

Already codified in `src/cyclerfinder/data/validate.py`'s module docstring
(#197): every `_LEVEL_EVIDENCE` gate is declared CONSISTENCY (same inputs,
different algorithm) or INDEPENDENCE (independently re-derived inputs); every
promotion requires at least one true independence gate; each cited gate carries
a one-line `shared with primary path:` declaration. Cross-linked here so the
whole doctrine is discoverable in one place.

## Per-interface external anchors (doctrine item 5)

The anchors this harness poisons ARE the per-interface external pins:

- epoch → `tests/verify/test_epoch_anchor.py` (JPL Horizons)
- frame / node convention → `tests/core/test_ephemeris_inclined.py` DE440 orbit-normal anchor
- sourced μ → `tests/core/test_cr3bp.py` (Ross & Roberts-Tsoukkas 2025)
