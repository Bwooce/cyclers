# `#855`: does `#839`'s C=3.13 heteroclinic connection promote `vaquero-31-c313-em-resonant-po-2013` to `cycler`?

**Task:** `#855`, registered 2026-08-21 (found during `#839`, not dispatched). Adjudicate whether
`#839`'s genuine, independently-verified `Wu(3:1) <-> Ws(2:1)` heteroclinic connection at C=3.13 —
touching `vaquero-31-c313-em-resonant-po-2013` via its OWN orbit as node31 — constitutes
"demonstrated transport utility" under schema v4.9/`#453`'s criterion, such that `orbit_class`
should move `resonant_po -> cycler`.

---

## Verdict (read this first)

**`orbit_class` stays `resonant_po`. No promotion, and this is not a close call once the
schema text and this project's own unbroken practice are read directly** — the honest residual
is that `resonant_po`'s *descriptive prose* now undersells this specific row, not that the
*enum value* is wrong. A comment-only `ADDED EVIDENCE` block was added to the row's
`orbit_class` field recording `#839`'s finding and this adjudication (Sec. 5); `orbit_class`
itself is unchanged.

## 1. What the schema actually says (read directly, not assumed)

`data/catalogue.schema.json`'s own v4.9 (`#453`) header text:

> adds a fifth `orbit_class` enum value `resonant_po` — a stable resonant / libration periodic
> orbit that is reachable at any epoch and repeats indefinitely (so it shares the cycler
> invariant `epoch_locked=false` / `n_returns='infinite'`) but has **NO demonstrated transport
> utility (it never encounters the secondary**; see Region B `#447`). Distinct from `cycler`,
> whose defining property IS transport between encounters.

`data/README.md`'s v4.9 entry is the same claim in fewer words: resonant_po is "carried for
known-class corroboration, not as a usable cycler — distinct from `cycler`, whose defining
property IS transport."

The subject of "it never encounters the secondary" is grammatically unambiguous: **the orbit**
(the antecedent of "it" is "a stable resonant/libration periodic orbit," not any object that
might later be found to connect to it). Nothing in either text mentions manifolds, transfers, or
connections. The criterion is about the periodic orbit's own trajectory.

## 2. This project's own practice, surveyed exhaustively (not a sample)

Every `orbit_class in {resonant_po, cycler}` row in `data/catalogue.yaml` was read (16 rows: 9
`resonant_po`, 7 `cycler`). Every single one determines `orbit_class` from exactly one
quantity — **the periodic orbit's own periselene** (its own closest approach to the Moon over
one period, from a dense propagation of the closed orbit itself) — compared against the lunar
SOI (66,182.9 km) and Hill radius (61,524.1 km):

| row | periselene (km) | vs SOI | class |
|---|---|---|---|
| `em-cycler-21-3d-spatial-2026` | 122,628 | 1.85x, outside | `resonant_po` |
| `casoliva-1-2c/1-2d/1-2e/2-1a/2-1b/3-2c` (6 rows) | 84,332–603,592 | 1.27x–9.12x, outside | `resonant_po` |
| `vaquero-21-c198-em-resonant-po-2013` | 86,911 | 1.31x, outside | `resonant_po` |
| **`vaquero-31-c313-em-resonant-po-2013` (this row)** | **66,995.2** | **1.012x, outside by 1.2%** | **`resonant_po`** |
| `casoliva-7-3a/7-3b/7-3c-em-cycler-2010` (3 rows) | 13,210–27,261 | 0.20x–0.41x, inside | `cycler` |
| `vaquero-21-c246/c247/c266-em-cycler-2013` (3 rows) | 17,675–49,700 | 0.27x–0.75x, inside | `cycler` |
| `vaquero-31-c254-em-cycler-2013` | 33,258 | 0.50x, inside | `cycler` |

Zero exceptions. No row's `orbit_class` comment cites a manifold leg, a heteroclinic connection,
or any object other than the orbit's own dense one-period propagation. `#811`'s own registration
for this exact 6-row Vaquero writeback states the discriminator explicitly: "`orbit_class`
determined PER MEMBER (not per family) against the lunar SOI... each cross-checked TWO
independent ways" — both of those two independent checks (`#799`'s DOP853 member_report and
`#811`'s own standalone Radau re-derivation) are propagations of the periodic orbit itself, never
of a transfer leg.

**This row's own `orbit_class` determination is itself an instance of the pattern**: `#811`'s
periselene value, 66,995.2 km (1.012x SOI, "OUTSIDE both, by only ~812 km"), is the sole basis
recorded on the row. `#839` never touched this quantity — it did not re-propagate the periodic
orbit, and its own note is explicit that the physical-character section is "reported as raw
physical data for `#855`'s adjudication... it is not itself a claim about `orbit_class`." The
periselene this row's classification rests on is therefore **unchanged** by `#839`'s work.

## 3. Direct precedent: this exact question was already answered once, on the sibling row

`vaquero-31-c254-em-cycler-2013` (the SAME 3:1 family's low-C end, already `cycler` because its
own periselene, 33,258 km, sits inside the SOI) carries `#828`'s `ADDED EVIDENCE` block recording
a verified heteroclinic connection (`#822`'s Wu(2:1)->Ws(3:1) free transfer) arriving at that
row. That block states, verbatim:

> Deliberately NOT a promotion: a heteroclinic connection is a two-object, same-model,
> same-fidelity TRANSPORT statement, whereas spec §14's ladder measures one object's own
> trajectory under increasing model fidelity.

That block is about `validation_level`, not `orbit_class` (c254 was already `cycler` on its own
periselene, so no `orbit_class` question arose there). But the underlying reasoning transfers
directly: a heteroclinic connection is evidence about a *two-object transport link*, categorically
distinct from evidence about *one object's own trajectory*. `#453`'s resonant_po/cycler criterion
is a one-object question ("it never encounters the secondary"). `#839`'s finding is two-object
evidence (a manifold leg connecting node31 to node21). Applying `#828`'s own distinction here
points the same way it did there: two-object transport evidence does not retroactively change
what is recorded about one object's own trajectory.

## 4. Corroborating evidence: the schema's own later design choices anticipated this exact case

Schema v5.3 (`#707`/`#708`, admitting the Uranus Umbriel-Titania CCR4BP torus-homoclinic
connection) explicitly considered and rejected reusing `resonant_po` for an object that DOES
have a verified connection but does NOT closely encounter the secondary moon on its own
trajectory:

> `orbit_class` enum widened with `torus_homoclinic` (chosen over reusing `quasi_cycler`, which
> presupposes a real multi-body encounter sequence this object doesn't have, and over
> `resonant_po`, **whose defining 'no demonstrated transport utility' property is the opposite of
> this row's whole point**)

Schema v5.4 (`#735`/`#736`, admitting the N=5 CRNBP torus) went further and declined even
`torus_homoclinic` for an object with a torus but no *computed* connection, precisely to avoid
overclaiming:

> `orbit_class` enum widened with `quasi_periodic_torus` (chosen over reusing `torus_homoclinic`,
> which presupposes a computed manifold connection this object doesn't have and would overclaim...
> and over `resonant_po`, whose strictly-periodic-orbit field conventions... cannot honestly
> encode a 2-D torus)

Both revisions demonstrate the project's own considered position, on the record twice: when an
object's evidentiary basis is a *connection* rather than the object's *own* orbital encounter
with the secondary, the correct response is a **new, honestly-scoped enum value**, not
reclassifying the object as `cycler` (which — per `#453` itself — asserts "transport between
encounters" as the orbit's OWN defining property) and not reusing `resonant_po` mislabelled
either. `vaquero-31-c313-em-resonant-po-2013` is not a new-class candidate on `#839`'s evidence
alone (one connection instance is not comparable in scope to the v5.3/v5.4 admissions, each of
which came with its own dedicated provenance block, multi-epoch real-ephemeris consistency work,
and a full design proposal) — but the schema's own precedent confirms that "verified connection,
no own-orbit encounter" is a recognized, previously-adjudicated combination, and this project's
answer to it has never been "call it `cycler`."

## 5. Independent re-verification of `#839`'s key numeric claim (not taken on say-so)

Per this project's standing discipline, `#839`'s selenocentric leg-minimum-radius numbers were
re-derived independently this task, not read out of `results.json` and trusted. Method: a
standalone inline DOP853 propagation (own `cr3bp_eom`, no `cyclerfinder` imports — same
independence discipline as `#811`'s own standalone Radau check), starting from `#839`'s own
recorded `crossing_state` for each of the two runs (`n_tau=48`, `n_tau=64`), integrated backward
for duration `t_u` and forward for duration `|t_s|` (both recorded in `results.json`), `rtol=atol
=1e-13`, 20,000-point dense sampling per leg, minimum selenocentric distance taken over the full
combined trajectory:

| run | `#839`'s value (km) | this task's independent re-derivation (km) | Δ (km) |
|---|---|---|---|
| `n_tau=48` | 46,247.64 | 46,247.57 | 0.06 |
| `n_tau=64` | 46,168.16 | 46,168.13 | 0.03 |

Agreement to <0.1 km on both — `#839`'s claim that the transfer leg dips to ~46,200 km
selenocentric (well inside the 66,182.9 km SOI, ~0.70x) is confirmed. **Scope of this
re-derivation, stated precisely**: it independently re-derives the minimum selenocentric radius
of the trajectory passing through `#839`'s own recorded `crossing_state` over `#839`'s own
recorded `t_u`/`t_s`. It inherits from `#839`'s own 11-test evidence battery the underlying claim
that this crossing state is a genuine, Newton-converged intersection of both manifolds (Newton
residual, ydot-sign gate, ghost guard, independent-Radau cross-check, forward/backward
re-approach — not re-run here, already independently re-verified by `#839`'s own two
phase-grid-resolution corroboration). This division matches the task's own framing: Sec. 1–4
above establish that this number, however solid, answers a two-object transport question, not
the one-object question `orbit_class` asks.

## 6. Conclusion

The tension the dispatch posed — periselene 66,995 km (1.2% outside SOI, `#811`, the orbit's own
trajectory) vs. transfer-leg minimum 46,200 km (well inside SOI, `#839`, a manifold leg reaching
the orbit) — resolves cleanly once the schema's actual subject ("it never encounters the
secondary," where "it" is the orbit) and this project's exhaustive, exception-free 16-row
practice are both consulted directly: `orbit_class` has always been, and remains, a
one-object question about the periodic orbit's own trajectory. `#839`'s finding is genuine,
independently corroborated (Sec. 5, and independently again by `#839`'s own two runs), and
scientifically important — but it is evidence about a different object (the transfer leg) than
the one `orbit_class` characterizes (the periodic orbit). `#811`'s own SOI-marginal periselene
determination for this row is untouched and still governs.

**Honesty note, stated rather than hidden**: `resonant_po`'s descriptive prose ("no demonstrated
transport utility") now genuinely undersells this specific row — a verified transport connection
exists at its own Jacobi constant, through its own orbit as node31. This is the same shape of
gap the row's own notes already flag for "stable" vs. this row's own UNSTABLE character
(|λ|=13.29): only `epoch_locked=false`/`n_returns='infinite'` is the schema-enforced invariant
(`tests/data/test_schema_v47_orbit_class.py`); the prose gloss is descriptive, not enforced, and
is imperfect here. If a future task wants to close that prose gap it should do so at the schema
level (e.g., an explicit "connection-adjacent" annotation convention, or refining `resonant_po`'s
description to acknowledge the boundary case) rather than by mislabelling this row `cycler` —
exactly the choice this project already made twice, in schema v5.3 and v5.4, when it faced the
same underlying tension on different objects.

## 7. Catalogue edit

Comment-only `ADDED EVIDENCE (#855 ...)` block appended to
`vaquero-31-c313-em-resonant-po-2013`'s `orbit_class` field comment (the field this evidence is
actually about, unlike `#828`/`#854`'s `validation_level`-field placement, which was correct for
their own validation-level question). `orbit_class` itself, `epoch_locked`, `n_returns`, and
`validation_level` are all unchanged — this is a comment-only edit, no semantic field value
changed.

## 8. Verification run

- Independent leg-minimum-radius re-derivation: Sec. 5, this task, ephemeral scratch script (not
  committed — verification only, no new artifact; the existing `#839` module and its 11-test
  battery already cover the underlying connection's correctness).
- Full `uv run pytest tests/data tests/search -q` ratchet (never a subset, per
  `[[feedback_catalogue_edits_run_all_ratchets]]`, `tests/search` split into two halves by
  filename per this environment's ~5-8 minute background-task ceiling): see commit log for exit
  status.
- `uv run ruff check .`, `uv run ruff format --check .`, full `uv run mypy src tests`: see commit
  log for exit status.
