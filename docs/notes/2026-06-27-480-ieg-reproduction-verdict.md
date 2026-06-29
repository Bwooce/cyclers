# #480 M1 Verdict — IEG EGGIE real-ephemeris reproduction (characterized NEGATIVE)

**Date:** 2026-06-27
**Task:** #480 Milestone 1 — faithful reproduction of the Hernandez-Jones-Jesick 2017
Io-Europa-Ganymede triple cycler (AAS 17-608) in real Galilean-moon ephemeris.
**Spec:** `docs/superpowers/specs/2026-06-26-480-resonance-lock-moontour-generator-design.md`
**Plan:** `docs/superpowers/plans/2026-06-27-480-ieg-real-eph-reproduction.md`
**Golden source (sourced invariants):** `docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md`
**Citation (verified-against-source, #486):** `anchor_for_key("hernandez-2017-ieg-608")`.

## Verdict in one line

The EGGIE 4-synodic ballistic cycler **does NOT reproduce** in real Galilean ephemeris
from the published-invariant seed. A correct Jupiter-central multiple-shooting corrector
**runs, makes genuine progress (defect cut ~5×), but converges into an off-paper,
high-energy, powered basin** — not the paper's 0.70 m/s ballistic family. This is a
**characterized negative gated on the publication gap** (the authors' unpublished state
vectors / V∞ directions), **not** a #473-style "the model is the wall" result.

## What was built (the M1 infrastructure — all landed, all green)

This milestone delivered the real-ephemeris Jovian substrate end-to-end, even though the
science verdict is negative:

| Unit | Deliverable | Commit |
|---|---|---|
| jup365 kernel helper | `verify/spice_kernels.py::ensure_jup365_kernel()` (path-or-RuntimeError, no auto-download) | `40d35cf` |
| Galilean real-eph backend | `core/ephemeris.py::_CentredSpiceBackend` (`model="spice"`, jup365 Io/Europa/Ganymede states, J2000-equatorial, center 599) | `c32b96e` |
| IEG seed adapter | `search/ieg_seed.py::ieg_eggie_seed()` — multi-rev Lambert legs, V∞-matched rev selection, `paper_departure_et()` | `924c416`, `7c4fead` |
| Real-eph correction (paper epoch) | Task 4 — IEG corrects against real Galilean ephemeris | `2eceba9` |
| Jupiter-central corrector | `nbody/jovian.py::jovian_shoot()` + `jovian_defect_residual()` (heliocentric `shooter.py` untouched, byte-for-byte) | `5c19006` |
| Verdict test + epoch scan | `tests/data/test_v4_jupiter.py`, `scripts/_ieg_rerun_scan_480.py` | `d5cf333`, `dbda56c` |

**Key architectural fact uncovered during Task 4:** `nbody/shooter.py::shoot()` is hardwired
**heliocentric** (MU_SUN + `PLANETS[body]` perturbers). The IEG tour is **Jupiter-central**.
Rather than touch the heliocentric shooter (and risk every existing cycler), Task 4a added a
parallel `jovian_shoot()` in `nbody/jovian.py` using `JovianRestrictedNBody` (MU_JUPITER,
jup365 moons). This is the reusable Jovian corrector for all future moon-tour work.

## The science result (numbers)

### Epoch scan (`scripts/_ieg_rerun_scan_480.py`, runlog `/tmp/ieg_rerun_scan_480.jsonl`)
- Scanned ±1 synodic period (±7.05 d) at 0.1 d, refined at 0.02 d, around the paper's
  2020-Oct-02 12:00 TDB departure (ET = 654912000.0).
- **Best `departure_et` = 654519744.0, offset −4.54 d** from the paper epoch. Sharp local
  minimum: **seed_defect_norm ≈ 3.89e5** there vs ~1.4e6 only 0.2 d away, and ~1.44e6 at
  the paper epoch itself. (The earlier single-rev scan's BEST_ET = 655335360.0, +4.9 d, is
  superseded — the correct multi-rev seed prefers −4.54 d.)

### Corrector verdict (`jovian_shoot`, multi-rev seed at BEST_ET)
- **CONVERGED: NO.** seed_defect_norm 3.89e5 → defect_norm **7.80e4** (ratio 0.20 — a
  genuine ~5× cut, real Gauss-Newton progress, not a stall at the seed).
- **Robust:** `max_nfev=120` stalls at nfev=34 at the *identical* fixed point as
  `max_nfev=30` (same defect, same V∞, same ΔV). The LM has reached its basin minimum;
  more iteration does not change the verdict.

| Encounter | Corrected V∞ (km/s) | Paper Table 4 (km/s) |
|---|---|---|
| Europa (depart) | 11.46 | 9.12 |
| Ganymede | 8.33 | 7.07 |
| Ganymede | 8.50 | 7.07 |
| Io | 4.89 | 8.38 |
| Europa (wrap) | 7.99 | 9.12 |

- **correction ΔV ≈ 5906 m/s** vs the paper's **0.70 m/s** ballistic close — four orders
  of magnitude apart. `bend_feasible = False` (corrected flyby bends not all achievable
  above the 100 km floor).

### Failure mode
The multi-rev Lambert seed, even at its best epoch, sits in a **different basin** (off-paper,
high-energy, powered) than the EGGIE ballistic cycler. The corrector descends to that
basin's minimum and stops; it never crosses to the paper's 0.70 m/s ballistic family.
Closing the gap needs a **family-targeted / better-topology seed** (the authors' actual
state vectors, or a homotopy from the ideal circular-coplanar EGGIE), **not** more
corrector iterations. The deep solve is FD-Jacobian-bound (~20 s/residual-eval) — the same
compute wall recorded for heliocentric multi-rev cyclers
(`memory/project_dsm_closure_modeljump_blocker.md`).

## Why this is the publication gap, not the #473 wall

#473 concluded "the model is the wall" for *single-window circular-coplanar* closes. This
result is different and more specific: the corrector here uses the **correct real-ephemeris
model** and a **multi-rev resonance-aware seed**, and it *does* find a converged fixed point
— just the wrong family. The barrier is **reconstruction underdetermination**: the paper
prints the sequence, leg ToFs, V∞ *magnitudes*, ΔV, and altitudes, but **not** the state
vectors or V∞ *directions*. The invariant set we can source is consistent with a continuum
of trajectories; ours lands in the high-energy one. This is exactly the **M1-stretch
publication-gap caveat** the spec called out as expected — invariant-match is the achievable
M1 bar, bit-exact author-trajectory match is the data-gated stretch.

## Catalogue / admission

**NO catalogue.yaml row was self-admitted.** There is no converged tour to admit — and even
a reproduced published tour would be a V4-ceiling known-reproduction (a published tour
cannot reach V5 by definition), reportable for separate **human** admission, never
self-admitted. Census/tier counts are unchanged.

## Tasks 5 / 6 disposition (reframed by the negative)

- **Task 5 (`v4_jupiter` gauntlet):** the V3/V4 gauntlet validates a *converged* tour. There
  is none, and the off-basin corrected state (5.9 km/s ΔV, bend-infeasible) is not the IEG
  cycler. Building gauntlet infra with no valid input and no positive control violates YAGNI
  and the positive-control discipline (`memory/feedback_verify_gauntlet_with_positive_control.md`).
  **Deferred** — built when a converged moon-tour exists (the natural consumer is the
  family-targeted-seed unblock below, and the pending moon-tour tasks #451 / #467 / #460).
- **Task 6 (invariant golden):** the convergence-gated reproduction assertion fails by
  design. The characterized negative is already affirmatively tested in
  `tests/data/test_v4_jupiter.py::test_ieg_jovian_shoot_at_best_epoch` (asserts *not*
  converged + defect improvement). The golden file `tests/verify/test_ieg_reproduction_golden.py`
  keeps the **convergence-independent citation gate active** and **skips** the reproduction
  assertion, pointing here. No tolerance was loosened to force a pass ("math decides").

## Follow-ups filed

1. **Family-targeted / homotopy seed for EGGIE** — the actual unblock: ramp ideal
   circular-coplanar EGGIE → real-eph, or acquire the authors' state vectors (M1-stretch,
   data-gated). This is what would cross from the off-paper basin to the 0.70 m/s family.
2. **`jovian_shoot` soft-wall-cap bug** — `max_wall_sec` is not enforced inside the scipy LM
   loop; probes overrun their budget. Cap inside the residual callback.
3. **`_CentredSpiceBackend` vs `nbody/jovian.py` Galilean `state()` duplication** — two
   jup365 readers now exist; reconcile to one provider.
4. **`v4_jupiter` gauntlet** — build when follow-up 1 yields a converged tour.
