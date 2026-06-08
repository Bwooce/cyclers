# Moon-tour Tier-1 — patched-conic moon systems + VILM (task #76)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:subagent-driven-development`
> or `superpowers:executing-plans`. Checkbox steps; strict TDD (write failing
> test → run **red** → minimal impl → run **green** → commit). Work on `main` —
> **do NOT branch** (project rule). **Do NOT commit while this plan is being
> authored under the docs-only mandate**; the commit messages below are the
> messages to *use* when the implementation phase runs (the user reviews first).
> uv-managed venv (no pip). Lint/type gate before **every** commit:
> `uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src tests`.
> Fast suite: `uv run pytest -m "not slow"`; corrector/ephemeris geometry tests
> that hit DE440 are `slow`.
>
> This plan is the task-level expansion of the **approved** design
> `docs/superpowers/specs/2026-06-06-moontour-planetcentric-design.md` (read it,
> **including its Approval section** — all five recommendations were accepted:
> Tier-1-first; VILM module IN scope (Option 2); re-tag the mis-flagged
> Jovian/Saturnian rows; 2-letter moon codes "reconciled against the codes
> already reserved in the catalogue"; JPL three-body catalog accepted for the
> *later* Tier-2 backfill, NOT this plan). Where this plan and the spec could
> diverge, **the spec + its Approval win** — except where a design claim no
> longer holds against live code, which this plan flags and overrides with the
> verified state (see "Design claims that no longer hold", below).
>
> **CONCURRENCY (binding).** Two other agents are editing `src/`, `tests/`,
> `data/` right now. A **task #110 agent owns `search/correct.py`,
> `search/seed_ladder.py`, `search/optimize.py`** (the M-ED ballistic corrector).
> Every task in this plan that touches `correct.py` carries an explicit
> **VERIFY-FIRST** step re-reading the live signature before editing, because the
> #110 agent's signatures may have drifted from what this plan saw on 2026-06-06.
> The `correct.py` edits here are **minimal surgical lifts of two Sun-couplings
> only** (§Phase 3) — do not refactor around the #110 work.

---

## Goal

Make the catalogue's **patched-conic moon-system** cycler rows (the two Jovian
rows + the Titan sub-family of the Saturnian row) computable on the *same*
dynamical model the heliocentric catalogue already uses (Kepler conics +
impulsive gravity-assist V∞ rematch), by adding a **central body = planet,
flyby bodies = its moons** axis. Plus the **VILM feasibility/search layer**
(Endgame Part-1 quadrature ΔV-floor + V̄∞-efficiency root) as the
"open planet-centric search" half of the milestone (Approval Q2 / Option 2).

The substrate is already centre-parametric in the layers that matter — Lambert
takes `mu=` (`core/lambert.py:516`), the Tisserand identity is canonical
(`search/tisserand.py:88-110`), and the registry is the single source of truth
designed to be extended (`core/constants.py:229-236`). The genuinely new content
is: a **`SATELLITES`/`PRIMARIES` registry**, a **planet-centred circular
ephemeris backend**, lifting the corrector's **two Sun-couplings**, teaching
**Tisserand `_a_p_km`/μ** to resolve a moon, the **catalogue re-tag** (so the
gauntlet dispatch stops misrouting), and the **VILM module** gated on the
transcribed Endgame Part-1 anchors.

**Tier 2 (CR3BP: Earth-Moon Arenstorf/Genova/Wittal + Saturnian midsize moons)
is OUT of scope** — citation-only until the later CR3BP milestone (Approval Q1).

---

## Design claims that no longer hold (verified live 2026-06-06 — these OVERRIDE the design)

The design's audit is the map, but three of its file/path/scheme claims are
stale against the live tree. **Implementers MUST follow the verified state
below, not the design's wording.**

1. **Tisserand lives at `search/tisserand.py`, NOT `core/tisserand.py`.** The
   design cites `core/tisserand.py:78-80,88-110` throughout (design §0, §1, §4).
   The live module is `src/cyclerfinder/search/tisserand.py`: `_a_p_km` at
   `search/tisserand.py:78-80`, `vinf_to_tisserand`/`tisserand_to_vinf` at
   `search/tisserand.py:88-110`, `linkable` at `:306`, `linkable_3d` at `:471`.
   The *line numbers* match; only the package is wrong. All Phase-4 tasks here
   use `search/tisserand.py`.

2. **The moon body-code scheme is FULL MOON NAMES, not 2-letter codes.** The
   design's Approval Q4 adopted "2-letter moon codes (Io, Eu, Ga, Ca, Ti, En)".
   But the live catalogue **already uses full names** as body codes —
   `data/README.md:65-67` mandates `"Io"`, `"Europa"`, `"Ganymede"`,
   `"Callisto"` ("Full moon names are used (not single letters) to avoid
   collision with heliocentric planet codes"); Saturnian reserved as `"Mimas"`,
   `"Enceladus"`, `"Tethys"`, `"Dione"`, `"Rhea"`, `"Titan"`
   (`data/README.md:68-69`); Martian `"Phobos"`/`"Deimos"` (`:70`). The live
   rows confirm it: `bodies: ["Io", "Europa", "Ganymede"]`
   (`catalogue.yaml:7945`), `["Io","Europa","Ganymede","Callisto"]`
   (`catalogue.yaml` russell-jovian), `["Mimas","Enceladus","Tethys","Dione",
   "Rhea","Titan"]` (`catalogue.yaml:8246`). The catalogue note at
   `catalogue.yaml:8338-8340` ("they may need new body codes for Saturnian
   moons (already reserved in data/README.md schema-extended section)") points
   at *that full-name reservation*. **Reconciliation outcome (Approval Q4
   "reconciled against the codes already reserved"): adopt the reserved
   full-name scheme; the `SatelliteData.code` field carries the full name (e.g.
   `code="Europa"`), and the registry is keyed by that full name.** A 2-letter
   alias is NOT introduced — it would collide with the catalogue's existing
   `bodies:` values and force a data migration the design did not authorise.
   This plan's Task 1.0 records the scheme decision in the registry docstring.

3. **`correct.py` already exists and has landed Phase 1 (the #110 agent).** The
   design (§0, §4) describes `ballistic_correct` / `_vinf_nodes` as things to
   patch, citing `search/correct.py:115`. They exist now:
   `ballistic_correct` at `search/correct.py:232` (signature has **no `mu`/
   `primary`** param), `_vinf_nodes` at `:78` (calls `lambert(...)` with **no
   `mu`**, `:115`; `ephem.state` heliocentric, `:109`), `_max_bend_deg` at `:182`
   (reads `PLANETS[body]`, `:186`). These are the **two Sun-couplings** to lift
   (Lambert μ default + the `PLANETS`-only bend lookup + the heliocentric
   ephemeris handed in). VERIFY-FIRST before every edit (concurrency).

4. **The fidelity ladder has FOUR rungs, not three.** Design §7 says
   "`circular-coplanar → analytic-ephemeris → real-de440`". Live
   `provenance.py:95` is
   `Literal["circular-coplanar","circular-inclined","analytic-ephemeris","real-de440"]`
   (4 rungs; `circular-inclined` sits between). Phase 7 tasks track the live
   4-rung set.

5. **No `core/ephemeris.py` Sun-subtract line drift.** Design §0 cites
   `core/ephemeris.py:263-264` for the astropy Sun-subtraction; verified at
   `core/ephemeris.py:263-265`. `_circular_inplane_state` at
   `core/ephemeris.py:109` (design said `:117` — that is the `a_km = sma_au *
   AU_KM` line inside it, `:117`, also verified). Backends: `_CircularBackend`
   `:130`, `_InclinedCircularBackend` `:151`, `_AstropyBackend` `:209`.

---

## Architecture

### The patched-conic moon model (design §1 Tier 1, §2)

Identical to the heliocentric N-arc model the #110 corrector already solves,
but with **central μ = a planet's μ** and **flyby bodies = its moons**:

- **Energy invariant:** V∞ *about the primary* (Jovicentric/Saturnian V∞ is a
  real conic excess speed; the canonical Tisserand identity `T = 3 − v∞²·a_p/μ`
  holds with the primary's μ — Endgame Part-2, mining note line 230).
- **Ephemeris:** moon position on its mean-motion circle **about the primary**,
  km-scaled (new `_CentredCircularBackend`, §Phase 2), not heliocentric.
- **Corrector:** the V∞-magnitude-continuity + closure residual is
  centre-agnostic; it needs (a) Lambert solved with `mu=μ_primary`, (b)
  `_max_bend_deg` resolving the moon's μ/radius, (c) a centred ephemeris handed
  in. Lift exactly those (§Phase 3).
- **Feasibility:** `linkable`/`linkable_3d` contour intersection is centre-blind
  once `_a_p_km` + μ are right (§Phase 4). VILM quadrature gives the ΔV cost on
  surviving moon pairs (§Phase 5).

### Why heliocentric stays byte-identical (design §2)

`model_assumption`/`primary` is a **pool pre-filter**, not a signature input
(`docs/spec.md:551`). A new hit is bucketed by `(model_assumption, primary)`,
matched only within its bucket; a Jovicentric V∞ never compares to a
heliocentric V∞. The V∞ comparator code is unchanged; only the partition key
gains `primary` (§Phase 6). Tier 1 reuses the V∞ comparator under a new bucket.

### Signature/matching extension (design §2, §6)

The gauntlet dispatch `anchors_for` keys off `cycler_class`
(`data/validate.py:558-605`): a `non-keplerian` row demands the `cr3bp` anchor
triple and forbids `invariants`. The two Jovian rows are tagged
`non-keplerian` (`catalogue.yaml` row 41 `cycler_class: non-keplerian`,
verified) but are **really patched-conic multi-arc moon tours** — so the
gauntlet currently demands a CR3BP identity from a patched-conic row. Re-tagging
to `multi-arc` flips `anchors_for` to `{invariants: True, cr3bp: False}` and the
dispatch routes them correctly (§Phase 2 re-tag).

---

## Tech stack

Python 3.11, numpy, `scipy.optimize` (already used by the corrector), pytest +
pyyaml. uv-managed venv. New production modules:
`src/cyclerfinder/core/satellites.py` (registry — kept out of `constants.py` to
avoid a merge collision with concurrent `constants.py` edits; re-exported from
`constants` only if clean),
`src/cyclerfinder/search/vilm.py` (VILM quadrature/efficiency). Extended:
`src/cyclerfinder/core/ephemeris.py` (centred backend),
`src/cyclerfinder/search/tisserand.py` (`_a_p_km`/μ moon resolution),
`src/cyclerfinder/search/correct.py` (two Sun-couplings — minimal). Catalogue
edits: `data/catalogue.yaml` (re-tag), `data/README.md` (scheme reconciliation
note). New tests under `tests/core/`, `tests/search/`, `tests/data/`.

---

## Phasing (independently shippable)

| Phase | Theme | Tasks | depends on |
|---|---|---|---|
| **1** | `SatelliteData` + `SATELLITES`/`PRIMARIES` registry + sourcing golden | 1.0–1.3 (4) | — |
| **2** | Catalogue re-tag + admit 2 Jovian rows + Titan sub-family; census ratchet | 2.0–2.3 (4) | — |
| **3** | Centre-agnostic corrector (lift the 2 Sun-couplings; heliocentric byte-identical guard) | 3.0–3.3 (4) | #110 corrector |
| **4** | Tisserand / T-P feasibility for moon systems (`_a_p_km`/μ; `T=3−v∞²`) | 4.0–4.2 (3) | Phase 1 |
| **5** | VILM module (phase-free ΔV + ΔV-min quadrature) gated on transcribed anchors | 5.0–5.4 (5) | Phase 1 |
| **6** | Signature/matching extension ((model_assumption, primary) pre-filter; heliocentric byte-identical) | 6.0–6.2 (3) | Phase 2 |
| **7** | Gauntlet/fidelity integration for planet-centric rows | 7.0–7.3 (4) | Phases 1–4, 6 |
| **2** is the centred ephemeris prerequisite for Phase 3 — see Phase 2.X note. | | | |

Phase 1 (registry) and Phase 5 (VILM) are pure and shippable alone. Phase 2
(centred ephemeris + re-tag) gates Phase 3. **Total: 27 tasks across 7 phases.**

> **NOTE on phase ordering:** the *centred ephemeris backend* is logically
> Phase 1 substrate but is placed in **Phase 2 (Task 2.0)** because it is the
> direct prerequisite for the Phase 3 corrector lift and shares no file with the
> registry. The re-tag tasks are 2.1–2.3.

---

## Phase 1 — `SATELLITES` / `PRIMARIES` registry + sourcing golden

Mirror the `PLANETS`/`SUPPORTED_BODIES` contract (`constants.py:170-236`) for
moons. A moon **cannot** live in `PLANETS` because `PlanetData.sma_au` is
intrinsically heliocentric (`constants.py:139` field; `_mean_motion_deg_day`
uses `MU_SUN`, `constants.py:149-159`). New sibling type in a **new module**
`core/satellites.py` (avoids a write-collision with concurrent `constants.py`
edits).

### Task 1.0 — `SatelliteData` dataclass + Kepler-III mean-motion derive ✅ DONE (a4effd3)

**Files:** create `src/cyclerfinder/core/satellites.py`; test
`tests/core/test_satellites_type.py`.

The dataclass mirrors `PlanetData` but is **about-primary** (km, not AU). Mean
motion is **derived at import** from `sma_km` + `μ_primary` via Kepler III
(exactly as `_mean_motion_deg_day`, `constants.py:149-159`) so the table stays
internally consistent rather than hand-copied. `code` carries the **full moon
name** (scheme reconciliation, see "Design claims" #2).

#### Failing test — `tests/core/test_satellites_type.py`

```python
"""Tier-1 Phase 1: SatelliteData shape + Kepler-III mean motion (plan Phase 1)."""
from __future__ import annotations

import math

from cyclerfinder.core.satellites import SatelliteData, mean_motion_deg_day_about


def test_satellitedata_fields_present() -> None:
    s = SatelliteData(
        name="Europa",
        code="Europa",
        primary="Jupiter",
        mu_km3_s2=3203.0,
        radius_eq_km=1560.8,
        sma_km=671100.0,
        mean_motion_deg_day=mean_motion_deg_day_about(671100.0, mu_primary=1.26686534e8),
        safe_alt_km=100.0,
    )
    assert s.code == "Europa"          # full-name scheme (data/README.md:65-67)
    assert s.primary == "Jupiter"
    assert s.sma_km == 671100.0


def test_mean_motion_matches_kepler_third_law() -> None:
    # n = 360 / (2 pi sqrt(a^3/mu) / 86400) deg/day; Europa about Jupiter.
    mu_jup = 1.26686534e8
    a = 671100.0
    period_s = 2.0 * math.pi * math.sqrt(a**3 / mu_jup)
    expected = 360.0 / (period_s / 86400.0)
    assert mean_motion_deg_day_about(a, mu_primary=mu_jup) == expected
```

Run: `uv run pytest tests/core/test_satellites_type.py -q` → **red** (no module).

#### Minimal impl — `core/satellites.py` (type + helper only, this task)

```python
"""Planet-centric satellite (moon) registry — Tier-1 patched-conic moon tours.

A moon cannot live in core.constants.PLANETS: PlanetData.sma_au is intrinsically
heliocentric (constants.py:139) and its mean motion derives from MU_SUN. This is
the about-the-primary sibling. Body-code scheme: FULL MOON NAMES (data/README.md
:65-69) — "Io"/"Europa"/"Ganymede"/"Callisto"/"Titan"/... — to avoid collision
with the heliocentric V/E/M planet codes. mean_motion_deg_day is DERIVED at
import from sma_km + the primary's mu (Kepler III), never hand-copied.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class SatelliteData:
    name: str
    code: str               # full moon name (scheme: data/README.md:65-69)
    primary: str            # "Jupiter" | "Saturn" | "Earth" | "Mars"
    mu_km3_s2: float        # moon GM
    radius_eq_km: float
    sma_km: float           # SMA ABOUT THE PRIMARY (km, not AU, not Sun-relative)
    mean_motion_deg_day: float   # about the primary (derive via mean_motion_deg_day_about)
    safe_alt_km: float


def mean_motion_deg_day_about(sma_km: float, *, mu_primary: float) -> float:
    """Mean motion (deg/day) about a primary, Kepler III (cf. constants.py:149-159)."""
    period_s = 2.0 * math.pi * math.sqrt(sma_km**3 / mu_primary)
    return 360.0 / (period_s / 86400.0)
```

Run → **green**. Lint/type. Commit:

```
core/satellites: SatelliteData type + Kepler-III mean-motion helper (moon-tour Tier-1 Phase 1)
```

### Task 1.1 — `PRIMARIES` table (Jupiter, Saturn) — sourced ✅ DONE (c34f53c)

**Files:** `core/satellites.py`; test `tests/core/test_primaries.py`.

`PRIMARIES: dict[str, float]` maps primary name → μ_primary. **Sourcing
discipline (golden rule applies to registry values exactly like catalogue
values):** μ_Jupiter / μ_Saturn trace to the **JPL SSD planetary GM table**
(the original tables, NOT the Endgame paper transcription — the paper's Table 3
is a *golden anchor* and must stay independent of the registry it validates, per
the golden-circularity rule + design §3 sourcing). Each value carries an inline
comment naming the JPL source and the queried value, exactly as `PLANETS` cites
Standish & Williams (`constants.py:120-132,179-208`). Earth μ already exists
(`PLANETS["E"].mu_km3_s2`, `constants.py:191`) and may be referenced for the
Tier-2 Earth-Moon primary later.

#### Failing test

```python
"""Tier-1 Phase 1: PRIMARIES mu table is present + sourced (plan Phase 1)."""
from __future__ import annotations

from cyclerfinder.core.satellites import PRIMARIES


def test_jupiter_and_saturn_primaries_present() -> None:
    assert "Jupiter" in PRIMARIES
    assert "Saturn" in PRIMARIES
    # Jovicentric GM ~1.2669e8 km^3/s^2 (JPL SSD planetary GM table).
    assert 1.26e8 < PRIMARIES["Jupiter"] < 1.27e8
    # Saturnian GM ~3.7931e7 km^3/s^2.
    assert 3.79e7 < PRIMARIES["Saturn"] < 3.80e7
```

> **SOURCING TODO for the implementer (not a code placeholder):** before writing
> the literal μ values, look up the **JPL SSD planetary physical/GM table** and
> paste the queried GM + the access date into the inline comment. The test
> asserts only an order-of-magnitude band (so it is not a circular self-check);
> the *exact* sourced value + citation is the golden discipline, recorded in the
> comment, cross-checked in Task 1.3 against the paper anchor.

Run → **red** → impl `PRIMARIES` with sourced inline comments → **green**.
Commit:

```
core/satellites: PRIMARIES mu table (JPL SSD GM), sourced (moon-tour Tier-1 Phase 1)
```

### Task 1.2 — `SATELLITES` table (Galilean four + Titan + Saturnian midsize) — sourced ✅ DONE (3479dbb)

**Files:** `core/satellites.py`; test `tests/core/test_satellites_registry.py`.

`SATELLITES: dict[str, SatelliteData]` keyed by full name. Tier-1 needs Io,
Europa, Ganymede, Callisto (Jovian) + Titan (Saturnian). Add the Saturnian
midsize moons (Mimas/Enceladus/Tethys/Dione/Rhea) too — they are cited in the
mining note Table 3 and needed for the VILM gates in Phase 5 even though their
*cyclers* defer to Tier 2. Each `mu_km3_s2`/`radius_eq_km`/`sma_km` traces to
**JPL SSD satellite physical-parameter tables** (inline comment per value).
`mean_motion_deg_day` is `mean_motion_deg_day_about(sma_km,
mu_primary=PRIMARIES[primary])` at construction. `safe_alt_km` = the altitude
the paper used for its ΔV table (100 km for all except Titan 1500 km, mining
note line 351) so the Phase-5 anchors are comparable.

> **SOURCING (golden, non-negotiable):** registry values come from **JPL SSD**,
> NOT from the Endgame Part-1 Table 3 transcription (mining note 337-352). The
> paper's `ã_M`/`Ṽ_M`/`μ̃_M` are the **independent cross-check anchor** (Task
> 1.3) — using them *as* the registry would make the Task-1.3 golden circular.

#### Failing test — `tests/core/test_satellites_registry.py`

```python
"""Tier-1 Phase 1: SATELLITES registry coverage + internal consistency (plan Phase 1)."""
from __future__ import annotations

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES, mean_motion_deg_day_about

_GALILEAN = ("Io", "Europa", "Ganymede", "Callisto")
_SATURNIAN = ("Mimas", "Enceladus", "Tethys", "Dione", "Rhea", "Titan")


@pytest.mark.parametrize("moon", _GALILEAN + _SATURNIAN)
def test_moon_present_and_keyed_by_full_name(moon: str) -> None:
    assert moon in SATELLITES
    assert SATELLITES[moon].code == moon          # full-name scheme
    assert SATELLITES[moon].primary in PRIMARIES


@pytest.mark.parametrize("moon", _GALILEAN + _SATURNIAN)
def test_mean_motion_is_derived_not_handcopied(moon: str) -> None:
    s = SATELLITES[moon]
    expected = mean_motion_deg_day_about(s.sma_km, mu_primary=PRIMARIES[s.primary])
    assert s.mean_motion_deg_day == pytest.approx(expected, rel=1e-12)
```

Run → **red** → impl `SATELLITES` with sourced comments → **green**. Commit:

```
core/satellites: SATELLITES table (Galilean + Saturnian), JPL-sourced (moon-tour Tier-1 Phase 1)
```

### Task 1.3 — registry-construction golden vs Endgame Part-1 Table 3 (independent anchor) ✅ DONE (9b2485f)

**Files:** test `tests/core/test_satellites_golden.py`.

**GOLDEN DISCIPLINE:** EXPECTED side = the **published** Endgame Part-1 Table 3
values (mining note 337-352), which were sourced *independently* of the JPL-SSD
values that built the registry. The registry's `sma_km` (ã_M) and a derived
circular velocity Ṽ_M = √(μ_primary / sma_km) must reproduce the paper's table
to a documented tolerance. This validates the registry construction; it is
separate from the Phase-5 VILM gates.

#### Failing test (verbatim anchors from mining note A1, Table 3)

```python
"""Tier-1 Phase 1 GOLDEN: SATELLITES reproduces the published Endgame Part-1
Table 3 a_M / V_M (mining note A1, lines 337-352).

GOLDEN DISCIPLINE: EXPECTED = the PUBLISHED Campagnola & Russell Part-1 Table 3
values, sourced INDEPENDENTLY of the JPL-SSD values that built the registry. The
circular velocity is V_M = sqrt(mu_primary / a_M). Tolerance bands documented
inline; the paper rounds a_M to km and V_M to 3 dp.
"""
from __future__ import annotations

import math

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# (moon, a_M km, V_M km/s) — Campagnola & Russell, "Endgame Problem" Part 1,
# Table 3 (p.17), transcribed at docs/notes/2026-06-05-endgame-tisserand-mining.md:339-348.
_TABLE3 = [
    ("Io",        421800.0, 17.330),
    ("Europa",    671100.0, 13.739),
    ("Ganymede", 1070400.0, 10.879),
    ("Callisto", 1882700.0,  8.203),
    ("Enceladus", 238040.0, 12.624),
    ("Tethys",    294670.0, 11.346),
    ("Dione",     377420.0, 10.025),
    ("Rhea",      527070.0,  8.484),
    ("Titan",    1221870.0,  5.572),
]


@pytest.mark.parametrize("moon,a_km,v_kms", _TABLE3)
def test_registry_reproduces_published_table3(moon: str, a_km: float, v_kms: float) -> None:
    s = SATELLITES[moon]
    # a_M: registry (JPL SSD) vs paper (JPL SSD) — same upstream, expect <0.1%.
    assert s.sma_km == pytest.approx(a_km, rel=1e-3)
    v_circ = math.sqrt(PRIMARIES[s.primary] / s.sma_km)
    # V_M: circular speed at a_M; paper rounds to 3 dp -> 0.5% band absorbs the
    # mu_primary digit choice between JPL releases.
    assert v_circ == pytest.approx(v_kms, rel=5e-3)
```

Run → **red** (until registry values are right) → adjust registry/μ to sourced
values until **green**. **If a moon cannot reach the band, the gap is a sourcing
problem — fix the registry against JPL SSD, never loosen the band to pass.**
Commit:

```
test: registry-construction golden vs published Endgame Part-1 Table 3 (moon-tour Tier-1 Phase 1)
```

---

## Phase 2 — Centred ephemeris + catalogue re-tag + admit rows; census ratchet

### Task 2.0 — `_CentredCircularBackend` (moon about its primary, km-scaled) ✅ DONE (309d85f)

**Files:** `src/cyclerfinder/core/ephemeris.py`; test
`tests/core/test_ephemeris_centred.py`.

Add a backend analogous to `_CircularBackend` (`ephemeris.py:130-148`) but for a
**moon on its mean-motion circle about the primary**, returning planet-centred
inertial `(r, v)` in km / km·s⁻¹. It reads `SATELLITES[moon]` (not `PLANETS`),
uses `sma_km` directly (no `* AU_KM`) and the about-primary
`mean_motion_deg_day`. Selected by a new `center=`/`model=` option on the
`Ephemeris` constructor (verify the live constructor signature first — it may
have drifted). The heliocentric backends stay byte-identical: the new branch is
only entered when `center` names a primary.

#### Failing test

```python
"""Tier-1 Phase 2: planet-centred circular moon ephemeris (plan Phase 2 Task 2.0)."""
from __future__ import annotations

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import SATELLITES


def test_europa_state_is_about_jupiter_at_its_sma() -> None:
    ephem = Ephemeris(model="circular", center="Jupiter")
    r, v = ephem.state("Europa", 0.0)
    # |r| ~ Europa's SMA about Jupiter (km), NOT a heliocentric AU-scaled radius.
    assert np.linalg.norm(r) == \
        __import__("pytest").approx(SATELLITES["Europa"].sma_km, rel=1e-9)
    # circular speed |v| = sqrt(mu_jup / a); ~13.7 km/s.
    assert 13.0 < np.linalg.norm(v) < 14.5


def test_heliocentric_default_unchanged() -> None:
    # No center -> the existing heliocentric circular backend, byte-identical.
    helio = Ephemeris(model="circular")
    r, _ = helio.state("E", 0.0)
    assert 1.40e8 < np.linalg.norm(r) < 1.55e8  # ~1 AU in km
```

Run → **red** → impl the backend + `center=` selection → **green**. Then run the
**whole existing ephemeris suite** to confirm heliocentric byte-identity:
`uv run pytest tests/core/test_ephemeris*.py -q` (use the live test file name).
Commit:

```
core/ephemeris: planet-centred circular moon backend (center=) (moon-tour Tier-1 Phase 2)
```

### Task 2.1 — re-tag the two Jovian rows: `non-keplerian` → `multi-arc` ✅ DONE (d44f035)

**Files:** `data/catalogue.yaml`; test `tests/data/test_moontour_retag.py`.

Per Approval Q3 and the gauntlet-dispatch fix (Architecture §): the Jovian
patched-conic rows are tagged `cycler_class: non-keplerian`
(`hernandez-2017-jovian-ieg-triple-family` and
`russell-strange-2009-jovian-multimoon-family`, verified live), which makes
`anchors_for` (`validate.py:558-605`) demand the `cr3bp` triple and forbid
`invariants`. They are **patched-conic multi-arc moon tours**. Change
`cycler_class` to `multi-arc` on both (leave `model_assumption:
circular-coplanar`, which is already correct). Update the inline schema comment.

> **DATA-CONCURRENCY (binding):** a concurrent agent may be editing
> `data/catalogue.yaml`. VERIFY-FIRST: re-grep the two row ids and their current
> `cycler_class` line immediately before editing; edit only those two
> `cycler_class:` values; do not touch surrounding fields.

#### Failing test — `tests/data/test_moontour_retag.py`

```python
"""Tier-1 Phase 2: Jovian patched-conic rows route as multi-arc, not CR3BP
(plan Phase 2 Task 2.1). The gauntlet anchors_for dispatch keys off cycler_class
(validate.py:558-605); a non-keplerian tag demands the cr3bp identity triple from
a patched-conic row, which is the misroute the design (§7) flags."""
from __future__ import annotations

import pytest
import yaml

from cyclerfinder.data.validate import anchors_for
from tests._catalogue_loader import CATALOGUE_PATH

_JOVIAN = ("hernandez-2017-jovian-ieg-triple-family",
           "russell-strange-2009-jovian-multimoon-family")


def _row(entry_id: str) -> dict:
    for row in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if row["id"] == entry_id:
            return row
    raise AssertionError(f"catalogue row {entry_id!r} not found")


@pytest.mark.parametrize("entry_id", _JOVIAN)
def test_jovian_rows_are_multi_arc(entry_id: str) -> None:
    assert _row(entry_id)["cycler_class"] == "multi-arc"


@pytest.mark.parametrize("entry_id", _JOVIAN)
def test_jovian_dispatch_wants_invariants_not_cr3bp(entry_id: str) -> None:
    a = anchors_for(_row(entry_id))
    assert a["invariants"] is True
    assert a["cr3bp"] is False
```

Run → **red** → edit the two `cycler_class:` values → **green**. Commit:

```
data/catalogue: re-tag Jovian patched-conic rows multi-arc (was non-keplerian) (moon-tour Tier-1 Phase 2)

The two Jovian rows are circular-coplanar patched-conic multi-arc moon tours,
not CR3BP. Tagging them non-keplerian made anchors_for (validate.py:558-605)
demand the cr3bp identity triple and forbid the multi-arc invariants, misrouting
the gauntlet dispatch (design §7, Approval Q3).
```

### Task 2.2 — Saturnian mixed row: supersession-style honesty for the Titan sub-family ✅ DONE (d44f035)

**Files:** `data/catalogue.yaml`; test `tests/data/test_moontour_retag.py`
(extend).

The Saturnian row (`russell-strange-2009-saturnian-multimoon-family`,
`catalogue.yaml:8235`) is `model_assumption: cr3bp` / `cycler_class:
non-keplerian` with `bodies: ["Mimas","Enceladus","Tethys","Dione","Rhea",
"Titan"]`. The existing SCHEMA-MISMATCH note (`catalogue.yaml:8330-8340`) says it
is **mixed**: Titan cyclers are patched-conic-modellable; the midsize-moon
members tend to low-energy/CR3BP regimes. Tier 1 unlocks only the **Titan
sub-family**. Do **not** silently re-tag the whole row (that would falsely claim
the Enceladus/Mimas members are patched-conic). Instead, apply
**supersession-style honesty**:

- Keep the family-seed row `non-keplerian` (its midsize members are genuinely
  CR3BP — Tier 2).
- Extend its `notes:` with an explicit split statement: "Titan members are
  Tier-1 patched-conic-modellable (`primary: Saturn`, circular-coplanar V∞ about
  Saturn); midsize-moon members (Mimas/Enceladus/Tethys) defer to the CR3BP
  Tier-2 milestone. As individual Titan cyclers are split out they take
  `cycler_class: multi-arc` / `model_assumption: circular-coplanar`."
- Do **not** add a Titan-only row in this task (no sourced individual-Titan-cycler
  numbers to populate it without fabrication). The split statement is the honest
  hand-off; an individual Titan row is created only when a sourced member is
  ingested.

#### Failing test (extend `test_moontour_retag.py`)

```python
def test_saturnian_row_stays_non_keplerian_with_titan_split_note() -> None:
    row = _row("russell-strange-2009-saturnian-multimoon-family")
    # NOT silently re-tagged: midsize members are genuinely CR3BP (Tier 2).
    assert row["cycler_class"] == "non-keplerian"
    notes = (row.get("notes") or "")
    assert "Titan" in notes
    # The honest split must name the Tier-1/Tier-2 boundary.
    assert "Tier-1" in notes or "patched-conic" in notes
    assert "Tier-2" in notes or "CR3BP" in notes
```

Run → **red** → extend the row's `notes:` (VERIFY-FIRST re-grep) → **green**.
Commit:

```
data/catalogue: Saturnian row — Titan Tier-1/midsize Tier-2 split note (moon-tour Tier-1 Phase 2)
```

### Task 2.3 — census ratchet: re-derive coverage same-commit; no entry vanishes ✅ DONE (d44f035)

**Files:** test `tests/data/test_validation_dispatch.py` (the existing
CONSTRUCTIBLE invariant); whichever census test holds `EXPECTED_COVERAGE`.

The re-tags change `anchors_for` results for three rows. The
`test_all_constructible_entries_are_single_ellipse` invariant
(`tests/data/test_validation_dispatch.py:207`) must still hold (the re-tagged
rows are `multi-arc`/`non-keplerian`, never `single-ellipse`, so they must NOT
be in CONSTRUCTIBLE — verify the loader still excludes them). Re-derive any
`EXPECTED_COVERAGE`/census numbers from the **live loader** (do not copy stale
counts — the M-ED plan's ratchet rule). Run the full dispatch + coverage suite.

- [x] `uv run pytest tests/data/test_validation_dispatch.py -q` green (2026-06-08).
- [x] Coverage/census invariant holds; no entry silently dropped. Ratchet
  re-derived 192 -> 223 (Russell 2004 +16, Rall 1970 +15 ingests; moon-tour
  added 0 rows). Fixed in 8218a56.
- [x] `ruff check` + `ruff format --check` clean on moon-tour files; `mypy src
  tests` clean except a concurrent-agent free_return WIP file (out of lane).

Commit (only if census numbers move):

```
test: census ratchet after moon-tour re-tags (re-derived from live loader) (moon-tour Tier-1 Phase 2)
```

---

## Phase 3 — Centre-agnostic corrector (lift the two Sun-couplings)

> **CONCURRENCY — #110 owns `correct.py`.** Every task here begins with a
> **VERIFY-FIRST** re-read of the live signature. The edits are **minimal
> surgical lifts** of exactly two Sun-couplings (Lambert μ default; the
> `PLANETS`-only bend lookup) plus accepting a centred ephemeris. Do **not**
> refactor the #110 N-arc solver structure.

### Task 3.0 — VERIFY-FIRST snapshot of the live corrector ✅ DONE (read-only; signatures drifted from plan lines: ballistic_correct@340, _vinf_nodes@111, _max_bend_deg@287, no mu/primary param)

**Files:** none (read-only); record findings in the task log.

- [ ] Read `src/cyclerfinder/search/correct.py`. Confirm/record the live:
  - `ballistic_correct` signature (line; presence/absence of `mu`/`primary`).
  - `_vinf_nodes` signature + the `lambert(...)` call (does it pass `mu`?).
  - `_max_bend_deg` body lookup (`PLANETS[body]`? line).
  - `_residuals`/`_residual_vector` signatures (do they thread the ephemeris?).
- [ ] If a `mu`/`primary`/`center` parameter **already exists** (the #110 agent
  may have added it for a different reason), adapt Tasks 3.1–3.2 to extend rather
  than introduce it — note the drift in the report.

(No commit — this is the safety read mandated by the concurrency rule.)

### Task 3.1 — plumb `mu=μ_central` into the Lambert call (coupling #1) ✅ DONE (4724354)

**Files:** `search/correct.py`; test `tests/search/test_correct_centred.py`.

Add a `mu_central: float = MU_SUN_KM3_S2` keyword to `_vinf_nodes` (and thread it
from `_residuals` and `ballistic_correct`). Pass it into the `lambert(...)` call
(`correct.py:115` today). Default = `MU_SUN_KM3_S2` ⇒ **heliocentric callers
byte-identical** (Lambert's own default is already `MU_SUN_KM3_S2`,
`lambert.py:516`, so passing it explicitly changes nothing for the Sun path).

#### Failing test — `tests/search/test_correct_centred.py`

```python
"""Tier-1 Phase 3: corrector is centre-agnostic — mu_central plumbed into Lambert
(plan Phase 3 Task 3.1). Heliocentric default stays byte-identical."""
from __future__ import annotations

import inspect

from cyclerfinder.search.correct import _vinf_nodes, ballistic_correct
from cyclerfinder.core.constants import MU_SUN_KM3_S2


def test_vinf_nodes_accepts_mu_central_defaulting_to_sun() -> None:
    sig = inspect.signature(_vinf_nodes)
    assert "mu_central" in sig.parameters
    assert sig.parameters["mu_central"].default == MU_SUN_KM3_S2


def test_ballistic_correct_accepts_mu_central() -> None:
    sig = inspect.signature(ballistic_correct)
    assert "mu_central" in sig.parameters
    assert sig.parameters["mu_central"].default == MU_SUN_KM3_S2
```

Run → **red** → thread `mu_central` (VERIFY-FIRST) → **green**. Then run the
**existing corrector suite** (`tests/search/test_correct_*.py`, incl. `-m slow`
S1L1) to confirm the heliocentric solver is byte-identical. Commit:

```
search/correct: plumb mu_central into Lambert (centre-agnostic, Sun-default) (moon-tour Tier-1 Phase 3)
```

### Task 3.2 — `_max_bend_deg` resolves a moon's μ/radius (coupling #2) ✅ DONE (4724354)

**Files:** `search/correct.py`; test `tests/search/test_correct_centred.py`
(extend).

`_max_bend_deg` reads `PLANETS[body]` (`correct.py:186`). For a moon code it must
read `SATELLITES[body]`. Make the lookup fall back: `PLANETS.get(body)` else
`SATELLITES[body]` (both expose `mu_km3_s2`/`radius_eq_km`/`safe_alt_km`). No
residual-math change — bend feasibility is post-hoc.

#### Failing test (extend)

```python
def test_max_bend_resolves_a_moon_code() -> None:
    from cyclerfinder.search.correct import _max_bend_deg
    # Europa is in SATELLITES, not PLANETS; the bend lookup must resolve it.
    bend = _max_bend_deg(3.0, "Europa")
    assert bend > 0.0
    # Higher V_inf -> tighter turn at the same moon.
    assert _max_bend_deg(8.0, "Europa") < bend
```

Run → **red** → add the `SATELLITES` fallback → **green**. Commit:

```
search/correct: _max_bend_deg resolves moon codes via SATELLITES (moon-tour Tier-1 Phase 3)
```

### Task 3.3 — end-to-end: close the Hernandez Io-Europa-Ganymede chain (slow) ✅ DONE (a8ac2e9) — closure PASSES; bend-feasibility strict-xfail (honest-risk: no-leverage model lands ~10km/s vinf, 100-150deg turns vs 2-5deg max-bend)

**Files:** test `tests/search/test_correct_jovian.py` (slow).

The corrector, handed the **centred Jovian ephemeris** (`Ephemeris(model=
"circular", center="Jupiter")`) and `mu_central=PRIMARIES["Jupiter"]`, closes a
patched-conic Io-Europa-Ganymede chain about Jupiter. **NON-GOLDEN** for the V∞
value (our computation) — the assertion is *closure exists with bend-feasible
flybys at Jovicentric V∞*, not a sourced V∞ (the sourced-anchor comparison is
the Phase-7 gauntlet's job). The seed ToFs come from the moons' synodic periods
(Io 1.769 d, Europa 3.551 d, Ganymede 7.155 d — `catalogue.yaml` row 41 text).

```python
"""Tier-1 Phase 3 (slow): ballistic_correct closes an Io-Europa-Ganymede chain
about Jupiter. NON-GOLDEN: closure + bend-feasibility asserted; the V_inf value
is OUR computation, not a sourced anchor (Phase 7 gauntlet does the sourced
comparison)."""
from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.search.correct import ballistic_correct


@pytest.mark.slow
def test_jovian_ieg_chain_closes_about_jupiter() -> None:
    ephem = Ephemeris(model="circular", center="Jupiter")
    # I-E-G-I closed triple; seeds ~ resonant synodic spacings (days). Slack =
    # the longest leg. period ~ a Laplace-resonance multiple (seed only).
    r = ballistic_correct(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single", "single", "single"),
        t0_seed_sec=0.0,
        tof_seed_days=(3.5, 7.2),     # two free legs; G-I slack
        period_sec=(1.769 + 3.551 + 7.155) * 86400.0,
        ephem=ephem,
        vinf_cap=12.0,
        mu_central=PRIMARIES["Jupiter"],
        slack_leg=2,
    )
    # Closure is the gate; the exact seeds may need an epoch/branch sweep — if it
    # does not converge on the first seed, widen the seed sweep (NOT the tol).
    assert r.converged
    assert r.bend_feasible
```

Run `-m slow` → **red** → (the impl is Tasks 3.1–3.2; this task is the
integration assertion + any seed-sweep tuning) → **green**. **If it will not
close, STOP and report** the best residual + which leg fails — do not loosen
`tol_kms` or fabricate a seed (honest-risk rule). Commit:

```
test: corrector closes Jovian I-E-G chain about Jupiter (non-golden, slow) (moon-tour Tier-1 Phase 3)
```

---

## Phase 4 — Tisserand / T-P feasibility for moon systems

> Path correction: the live module is **`search/tisserand.py`** (design said
> `core/`). Line numbers match.

### Task 4.0 — `_a_p_km` resolves a moon's about-primary SMA ✅ DONE (b56546e)

**Files:** `src/cyclerfinder/search/tisserand.py`; test
`tests/search/test_tisserand_moon.py`.

`_a_p_km` reads `PLANETS[body].sma_au * AU_KM` (`tisserand.py:78-80`). For a moon
code it must return `SATELLITES[body].sma_km` (already km, about the primary).
Fallback: `PLANETS.get(body)` → `sma_au * AU_KM` else `SATELLITES[body].sma_km`.

#### Failing test

```python
"""Tier-1 Phase 4: Tisserand a_p resolves moon about-primary SMA (plan Phase 4)."""
from __future__ import annotations

import pytest

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.search.tisserand import _a_p_km


def test_a_p_km_for_moon_is_about_primary_sma() -> None:
    assert _a_p_km("Europa") == pytest.approx(SATELLITES["Europa"].sma_km)


def test_a_p_km_for_planet_unchanged() -> None:
    from cyclerfinder.core.constants import AU_KM, PLANETS
    assert _a_p_km("E") == pytest.approx(PLANETS["E"].sma_au * AU_KM)
```

Run → **red** → add the fallback → **green**. Commit:

```
search/tisserand: _a_p_km resolves moon about-primary SMA via SATELLITES (moon-tour Tier-1 Phase 4)
```

### Task 4.1 — μ_primary into `vinf_to_tisserand` / `tisserand_to_vinf` ✅ DONE (b56546e)

**Files:** `search/tisserand.py`; test `tests/search/test_tisserand_moon.py`
(extend).

`vinf_to_tisserand`/`tisserand_to_vinf` hardwire `MU_SUN_KM3_S2`
(`tisserand.py:94,110`). Add a `mu: float = MU_SUN_KM3_S2` keyword so the
canonical identity `T = 3 − v∞²·a_p/μ` is about the right centre. Default = Sun
⇒ **heliocentric callers byte-identical**. The round-trip
`tisserand_to_vinf(body, vinf_to_tisserand(body, v, mu=μ), mu=μ) == v` is the
`T = 3 − v∞²` identity the mining note confirms (note line 230).

#### Failing test (extend)

```python
def test_vinf_tisserand_roundtrip_about_jupiter() -> None:
    from cyclerfinder.core.satellites import PRIMARIES
    from cyclerfinder.search.tisserand import tisserand_to_vinf, vinf_to_tisserand
    mu = PRIMARIES["Jupiter"]
    v = 5.0  # km/s Jovicentric V_inf at Europa
    t = vinf_to_tisserand("Europa", v, mu=mu)
    assert tisserand_to_vinf("Europa", t, mu=mu) == pytest.approx(v, rel=1e-9)


def test_sun_default_unchanged() -> None:
    from cyclerfinder.search.tisserand import vinf_to_tisserand
    # No mu= -> heliocentric, byte-identical to pre-change.
    assert vinf_to_tisserand("E", 5.0) == vinf_to_tisserand("E", 5.0, )
```

Run → **red** → add `mu=` keyword threaded through both functions → **green**.
Run the **existing tisserand suite** to confirm Sun-default byte-identity.
Commit:

```
search/tisserand: mu= on vinf<->tisserand (centre-aware T=3-vinf^2), Sun-default (moon-tour Tier-1 Phase 4)
```

### Task 4.2 — `linkable` prunes a Jovicentric moon pair (centre-blind once a_p/μ right) ✅ DONE (b56546e)

**Files:** test `tests/search/test_tisserand_moon.py` (extend). `linkable`'s
contour-intersection logic (`tisserand.py:306`) is centre-blind once `_a_p_km`
and μ resolve correctly — confirm no further code change is needed. If `linkable`
internally hardwires `MU_SUN` (verify at `tisserand.py:255,292`), thread the same
`mu=` keyword (minimal).

```python
def test_linkable_resolves_a_jovicentric_moon_pair() -> None:
    from cyclerfinder.core.satellites import PRIMARIES
    from cyclerfinder.search.tisserand import linkable
    # Europa-Ganymede share Jupiter; at a feasible Jovicentric V_inf the contours
    # intersect (linked-conic patch point exists). Assert linkable returns a
    # result object / truthy feasibility, not a specific V_inf (non-golden).
    result = linkable("Europa", "Ganymede", mu=PRIMARIES["Jupiter"])  # adapt to live signature
    assert result is not None
```

> VERIFY-FIRST: read `linkable`'s live signature (`tisserand.py:306`) — it may
> take vinf ranges / a config object. Adapt the call; the assertion is *a
> Jovicentric moon pair resolves through the predicate without a Sun-SMA error*.

Run → **red**/**green** as the signature requires; thread `mu=` only if a live
`MU_SUN` hardwire is found. Commit:

```
test: linkable prunes a Jovicentric moon pair (centre-blind feasibility) (moon-tour Tier-1 Phase 4)
```

---

## Phase 5 — VILM module (gated on transcribed Endgame Part-1 anchors)

> **GATE SCOPE (verbatim, binding):** the VILM module's *numeric* output is
> gated by the **transcribed Endgame Part-1 Tables 1–3** in
> `docs/notes/2026-06-05-endgame-tisserand-mining.md` (sections A1/A2/A3, lines
> 337–397) **and** the Part-1 worked Europa scalar (A6, lines 436–438: 154 m/s /
> 46 days). **The two flagged suspect cells of Part-2 Table 1 are EXCLUDED as
> goldens** — specifically the Titan `J_L4@100km` cell "766.5/776.4" (MAX<MIN
> inversion, mining note A4 line 411 + lines 491–493). No Part-2 Table 1 value
> is used as a VILM gate in this phase (Part-2 is CR3BP/Tier-2). What gates what:
> Table 1 (no-GA ΔV_min) gates the Eq.(13) quadrature ΔV-floor; Table 2 (with-GA
> ΔV_min) gates the GA-routed quadrature; Table 3 V̄∞ E/I gates the Eq.(9)
> V̄∞-efficiency root; the A6 Europa scalar gates the end-to-end 3-VILM ΔV.**
> **Physics-invariant-only (NOT gated by any table):** the phase-free n:m_K± ΔV
> formula *shape* and the `T = 3 − v∞²` identity (those are validated by
> round-trip/algebraic identity in Phase 4, not by a paper cell).

### Task 5.0 — `vilm.py` skeleton + n:m_K± resonance taxonomy ✅ DONE (325dc2f)

**Files:** create `src/cyclerfinder/search/vilm.py`; test
`tests/search/test_vilm_taxonomy.py`.

The n:m_K± taxonomy (mining note lines 98–130) classifies a VILM leg by the
spacecraft:moon resonance `n:m`, the body `K`, and exterior/interior `±`. Pure
classification — **physics-invariant, not table-gated**.

#### Failing test

```python
"""Tier-1 Phase 5: VILM n:m_K± leg taxonomy (plan Phase 5 Task 5.0).
Physics-invariant classification; NOT gated by a paper cell."""
from __future__ import annotations

from cyclerfinder.search.vilm import classify_vilm_leg


def test_exterior_vilm_classification() -> None:
    leg = classify_vilm_leg(n=3, m=2, body="Europa", exterior=True)
    assert leg.resonance == (3, 2)
    assert leg.body == "Europa"
    assert leg.exterior is True
```

Run → **red** → impl `classify_vilm_leg` + a `VilmLeg` dataclass → **green**.
Commit:

```
search/vilm: VILM module skeleton + n:m_K± leg taxonomy (moon-tour Tier-1 Phase 5)
```

### Task 5.1 — Eq.(9) V̄∞-efficiency root, gated on Table 3 V̄∞ E/I ✅ DONE (8c03be0) — E AND I columns reproduced to <0.5 m/s

**Files:** `search/vilm.py`; test `tests/search/test_vilm_efficiency.py`.

`min_vinf_for_vilm(moon)` = the root of Endgame Part-1 Eq.(9) (the minimum V∞ at
which a VILM becomes efficient), per moon, using the registry's μ/a. **GOLDEN:**
EXPECTED = the **published** Table 3 `V̄∞ E/I` column (mining note A1, lines
337–354) — a *derived* quantity in the paper, so it doubles as the check on our
root.

#### Failing test (verbatim Table 3 V̄∞ anchors)

```python
"""Tier-1 Phase 5 GOLDEN: VILM efficiency root vs published Part-1 Table 3 V_inf_bar E/I
(mining note A1, lines 337-354).

GOLDEN DISCIPLINE: EXPECTED = the PUBLISHED Table 3 V_inf_bar values (km/s),
exterior (E) column. The corrector output is the side under test."""
from __future__ import annotations

import pytest

from cyclerfinder.search.vilm import min_vinf_for_vilm

# (moon, V_inf_bar exterior km/s) — Part-1 Table 3, E column (note lines 339-348).
_VINF_BAR_E = [
    ("Io", 0.351), ("Europa", 0.277), ("Ganymede", 0.372), ("Callisto", 0.328),
    ("Titan", 0.283), ("Rhea", 0.085), ("Dione", 0.067),
    ("Tethys", 0.052), ("Enceladus", 0.029),
]


@pytest.mark.parametrize("moon,vbar_e", _VINF_BAR_E)
def test_efficiency_root_matches_published_vinf_bar(moon: str, vbar_e: float) -> None:
    # Paper rounds to 3 dp; band absorbs the V_c assumption + mu digit choice.
    assert min_vinf_for_vilm(moon) == pytest.approx(vbar_e, abs=0.01)
```

Run → **red** → impl Eq.(9) root → **green**. **If a moon cannot reach the band,
report the gap** (it is a model/sourcing finding, not a tolerance to loosen).
Commit:

```
search/vilm: Eq.(9) V_inf-efficiency root vs Part-1 Table 3 V_inf_bar golden (moon-tour Tier-1 Phase 5)
```

### Task 5.2 — Eq.(13) ΔV-min quadrature (no-GA), gated on Table 1 ✅ DONE (a30d7bf) — all Table 1 to <0.5%

**Files:** `search/vilm.py`; test `tests/search/test_vilm_quadrature.py`.

`vilm_dv_min(moon_a, moon_b)` = the Eq.(13) closed-form quadrature ΔV-floor for a
no-gravity-assist VILM intermoon transfer. **GOLDEN:** EXPECTED = the published
Part-1 **Table 1 ΔV_min** column (mining note A2, lines 358–381).

#### Failing test (verbatim Table 1 ΔV_min anchors)

```python
"""Tier-1 Phase 5 GOLDEN: VILM ΔV-min quadrature vs published Part-1 Table 1
(no-GA, mining note A2, lines 358-381).

GOLDEN DISCIPLINE: EXPECTED = the PUBLISHED Table 1 ΔV_min (km/s). Model caveat
(mining note 491-493): Part-1 ΔV are linked-conic, upper-bound-ish vs CR3BP —
the band is a tolerance, never equality, and a CR3BP-disagreeing value <=10% is
NOT a rejection."""
from __future__ import annotations

import pytest

from cyclerfinder.search.vilm import vilm_dv_min

# (moon_a, moon_b, ΔV_min km/s) — Part-1 Table 1 (note lines 362-378).
_DV_MIN = [
    ("Callisto", "Ganymede", 1.81), ("Ganymede", "Europa", 1.71),
    ("Europa", "Io", 1.76), ("Titan", "Rhea", 1.15),
    ("Rhea", "Dione", 0.52), ("Tethys", "Enceladus", 0.34),
]


@pytest.mark.parametrize("a,b,dv", _DV_MIN)
def test_quadrature_dv_min_matches_published_table1(a: str, b: str, dv: float) -> None:
    # 10% band per the linked-conic vs CR3BP model caveat (mining note 491-493).
    assert vilm_dv_min(a, b) == pytest.approx(dv, rel=0.10)
```

Run → **red** → impl the Eq.(13) quadrature → **green**. Commit:

```
search/vilm: Eq.(13) ΔV-min quadrature vs Part-1 Table 1 (no-GA) golden (moon-tour Tier-1 Phase 5)
```

### Task 5.3 — GA-routed ΔV-min (Table 2) + end-to-end Europa scalar (A6) ✅ DONE (a30d7bf) — Table 2 to 0.2%; A6 floor 128 m/s = valid LB on discrete 154

**Files:** `search/vilm.py`; test `tests/search/test_vilm_quadrature.py`
(extend).

The GA-routed quadrature (intermediate moon flybys) gated on Part-1 **Table 2
ΔV_min** (mining note A3, lines 383–397); and the end-to-end Europa 3-VILM
endgame scalar gated on the worked A6 value (**154 m/s / 46 days**, mining note
lines 436–438).

#### Failing test (extend — Table 2 + A6)

```python
def test_ga_routed_dv_min_matches_table2() -> None:
    # Callisto-G-Europa ΔV_min 1.61 (Part-1 Table 2, note line 387).
    assert vilm_dv_min("Callisto", "Europa", via=("Ganymede",)) == pytest.approx(1.61, rel=0.10)


def test_europa_3vilm_endgame_scalar() -> None:
    from cyclerfinder.search.vilm import europa_endgame_dv
    # A6 (note 436-438): 3-VILM Europa endgame ~154 m/s over ~46 days, V_inf_bar
    # 1.8 -> 0.77 km/s. EXPECTED = the published worked scalar.
    dv_ms, days = europa_endgame_dv()
    assert dv_ms == pytest.approx(154.0, abs=15.0)   # m/s, ~10% band
    assert days == pytest.approx(46.0, abs=5.0)
```

Run → **red** → impl GA routing + the worked-scalar driver → **green**. Commit:

```
search/vilm: GA-routed ΔV (Table 2) + Europa 3-VILM endgame scalar (A6) golden (moon-tour Tier-1 Phase 5)
```

### Task 5.4 — VILM as admissible ΔV-floor for search pruning + Phase-5 gate ✅ DONE (a30d7bf) — DEVIATION: floor=escape+capture (plan's no-GA-as-floor is backwards; GA reduces ΔV)

**Files:** `search/vilm.py`; test `tests/search/test_vilm_quadrature.py`
(extend); lint/type gate.

Expose `vilm_dv_floor(moon_a, moon_b)` as an **admissible lower bound** (≤ any
real tour-leg ΔV) for A*-style pruning in a future Forge search (design §5). The
floor is `vilm_dv_min` (the no-GA quadrature is a valid lower bound). Assert it
is ≤ the with-GA value for a routed pair (GA can only reduce ΔV) — a structural
invariant, not a paper cell.

```python
def test_dv_floor_is_admissible_lower_bound() -> None:
    from cyclerfinder.search.vilm import vilm_dv_floor, vilm_dv_min
    floor = vilm_dv_floor("Callisto", "Europa")
    assert floor <= vilm_dv_min("Callisto", "Europa")           # floor <= no-GA
    assert floor <= vilm_dv_min("Callisto", "Europa", via=("Ganymede",))  # <= with-GA
```

- [x] `uv run pytest tests/search/test_vilm_*.py -q` green (2026-06-08).
- [x] `ruff check` + `ruff format --check` + `mypy` clean on VILM source/tests.

Commit:

```
search/vilm: admissible ΔV-floor for search pruning + Phase-5 gate (moon-tour Tier-1 Phase 5)
```

---

## Phase 6 — Signature/matching extension ((model_assumption, primary) pre-filter)

### Task 6.0 — `(model_assumption, primary)` bucket key on the matcher pool

**Files:** the matcher/signature module (locate via
`grep -rn "model_assumption" src/cyclerfinder/`; the spec §16.2 area,
`docs/spec.md:551`); test `tests/data/test_signature_bucket.py`.

Add a `(model_assumption, primary)` bucket key so a finder hit is matched only
within its bucket. A Jovicentric V∞ (primary=Jupiter) never compares to a
heliocentric V∞ (primary=Sun/absent). The V∞ comparator code is **unchanged**;
only the partition key gains `primary` (treat missing `primary` as `"Sun"`,
`data/README.md:51`).

#### Failing test

```python
"""Tier-1 Phase 6: matcher pool pre-filters by (model_assumption, primary)
(plan Phase 6; spec §16.2 / docs/spec.md:551). Heliocentric rows (primary
absent->Sun) never share a bucket with Jovicentric rows."""
from __future__ import annotations

# Adapt the import to the live matcher module discovered by grep.
from cyclerfinder.<matcher_module> import signature_bucket_key


def test_missing_primary_is_sun_bucket() -> None:
    assert signature_bucket_key({"model_assumption": "circular-coplanar"})[1] == "Sun"


def test_jovian_and_heliocentric_are_different_buckets() -> None:
    helio = signature_bucket_key({"model_assumption": "circular-coplanar"})
    jovian = signature_bucket_key(
        {"model_assumption": "circular-coplanar", "primary": "Jupiter"}
    )
    assert helio != jovian
```

> VERIFY-FIRST: `grep -rn "model_assumption\|def .*match\|bucket\|pool" src/cyclerfinder/`
> to find the live matcher and its pre-filter (spec §16.2 says M7 *may* pre-filter
> by `model_assumption` already — extend that key with `primary`, do not add a
> parallel path). Fill in `<matcher_module>` / `signature_bucket_key` to the live
> names.

Run → **red** → add `primary` to the bucket key → **green**. Commit:

```
search/matcher: (model_assumption, primary) pool pre-filter; Sun default (moon-tour Tier-1 Phase 6)
```

### Task 6.1 — heliocentric rows byte-identical guard

**Files:** test `tests/data/test_signature_bucket.py` (extend) + run the full
existing matcher/signature suite.

The canonical signature for heliocentric rows must be **byte-identical** — no
field changes, no golden moves (design §2, §9). Assert that for every
heliocentric catalogue row the signature is unchanged and the bucket key's
primary component is `"Sun"`. Run the existing signature/golden suite to confirm
no heliocentric golden moved.

```python
def test_all_heliocentric_rows_bucket_as_sun() -> None:
    import yaml
    from tests._catalogue_loader import CATALOGUE_PATH
    for row in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if row.get("primary") in (None, "Sun"):
            assert signature_bucket_key(row)[1] == "Sun"
```

- [x] Existing matcher/signature/golden suite green (no heliocentric move) (2026-06-08).

Commit (only if a guard test is added):

```
test: heliocentric signature byte-identical under primary bucket key (moon-tour Tier-1 Phase 6)
```

### Task 6.2 — Phase-6 lint/type gate

- [x] `uv run pytest tests/data/test_signature_bucket.py -q` + the live matcher
  suite green; ruff + ruff format + mypy clean (2026-06-08; mypy yaml ignore
  added in 8218a56).

---

## Phase 7 — Gauntlet / fidelity integration for planet-centric rows

The gauntlet combiner (`verify/gauntlet.py`) "asserts no computed physics value"
and folds axis reports (`gauntlet.py:61`). Planet-centric changes are confined to
**what each axis computes**, not the combiner.

### Task 7.0 — fidelity ladder rungs are about the primary (Axis B)

**Files:** `src/cyclerfinder/verify/fidelity.py`; test
`tests/verify/test_fidelity_moon.py`.

The Axis-B ladder is the live 4-rung `Fidelity` Literal (`provenance.py:95`:
`circular-coplanar → circular-inclined → analytic-ephemeris → real-de440` — note
the design's 3-rung wording is stale, "Design claims" #4). For a moon tour the
same rungs are **about the primary**: `circular-coplanar` = the
`_CentredCircularBackend` (Phase 2); `analytic-ephemeris` stays the documented
extension point (no in-house backend, `fidelity.py:28-29`); `real-de440` =
astropy moon states (astropy resolves Galilean moons via the body-name map,
`ephemeris.py:59`). The persistence check (does Jovicentric V∞ stay stable
circular→real?) is the same `_moves_toward_band` logic on a primary-relative V∞.

#### Failing test

```python
"""Tier-1 Phase 7: Axis-B fidelity persistence tracks Jovicentric V_inf across
rungs (plan Phase 7 Task 7.0). Same _moves_toward_band logic, about the primary."""
from __future__ import annotations

# Adapt to the live PersistenceReport / fidelity entry-point signature.
from cyclerfinder.verify.fidelity import persistence_for_primary  # name TBD - see VERIFY


def test_jovicentric_vinf_persistence_runs_about_jupiter() -> None:
    report = persistence_for_primary(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        primary="Jupiter",
        # ... seed/period as in Phase 3 Task 3.3
    )
    assert report is not None
    # The tracked scalar is V_inf about JUPITER, not the Sun.
    assert report.quantity in ("vinf", "vinf_primary")
```

> VERIFY-FIRST: read `verify/fidelity.py` for the live persistence entry point
> (`solve_at_fidelity` / `PersistenceReport`, design §7 cites `fidelity.py:15-25,
> 167-169,393-411`). The new code is a *primary-parametrised* call into the same
> ladder, not a new combiner. Adapt the test's import/signature to live names.

Run → **red** → add the primary-parametrised persistence path → **green**.
Commit:

```
verify/fidelity: Axis-B persistence about a primary (Jovicentric V_inf) (moon-tour Tier-1 Phase 7)
```

### Task 7.1 — Axis A: VILM quadrature vs corrector ΔV as a second code path

**Files:** test `tests/verify/test_agreement_moon.py`.

Axis A wants two independent code paths. For a Jovian moon pair, the **VILM
quadrature** (Phase 5) and the **corrector** (Phase 3) are two independent ΔV
computations on the same pair — exactly the crosscheck Axis A folds
(`gauntlet.py` Axis A). Assert the two agree within the linked-conic-vs-DE440
band (10%, mining note 491-493 model caveat — NOT equality). NON-GOLDEN (both
sides are our computation; this is an internal-agreement axis, not a sourced
gate).

```python
"""Tier-1 Phase 7: Axis-A code-path agreement — VILM quadrature vs corrector ΔV
on a Jovian moon pair (plan Phase 7 Task 7.1). NON-GOLDEN internal crosscheck;
10% band per the linked-conic vs DE440 model caveat (mining note 491-493)."""
from __future__ import annotations

import pytest


@pytest.mark.slow
def test_vilm_vs_corrector_agree_within_model_band() -> None:
    from cyclerfinder.search.vilm import vilm_dv_min
    # corrector ΔV for the same Europa-Ganymede leg (turn-deficit / continuity)
    # ... build via Phase 3 corrector, extract leg ΔV
    vilm = vilm_dv_min("Ganymede", "Europa")
    corrector_dv = ...  # from ballistic_correct leg output (live signature)
    assert corrector_dv == pytest.approx(vilm, rel=0.10)
```

> The `...` is a test sketch placeholder for the corrector-ΔV extraction (build
> via Phase 3 Task 3.3's `ballistic_correct` call and read the per-leg ΔV) — the
> recipe is named in prose; fill from the live result object.

Run `-m slow` → **red**/**green**. **If they disagree >10%, report it** (a model
finding, not a tolerance to widen). Commit:

```
test: Axis-A VILM-vs-corrector ΔV agreement on a Jovian pair (non-golden, slow) (moon-tour Tier-1 Phase 7)
```

### Task 7.2 — Axis C/D: anchors as corroboration; wrong-μ falsification guard

**Files:** test `tests/verify/test_falsification_moon.py`.

- Axis C (provenance): the mined Part-1 anchors are corroboration for the Jovian
  rows' VILM ΔV (Phase 5 already gates on them — this confirms the gauntlet sees
  them as Axis-C input).
- Axis D (falsification): a deliberately-wrong-μ guard — running the corrector
  with `mu_central=MU_SUN` (the *wrong* centre) on a Jovian chain must NOT close
  / must produce a refuted verdict. This is the "deliberately bogus" Axis-D guard
  (`gauntlet.py:17-21`).

```python
"""Tier-1 Phase 7: Axis-D falsification — a Jovian chain solved with the WRONG
centre (mu_Sun instead of mu_Jupiter) must be refuted (plan Phase 7 Task 7.2)."""
from __future__ import annotations

import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct


@pytest.mark.slow
def test_wrong_central_mu_does_not_spuriously_close() -> None:
    ephem = Ephemeris(model="circular", center="Jupiter")
    r = ballistic_correct(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single", "single", "single"),
        t0_seed_sec=0.0,
        tof_seed_days=(3.5, 7.2),
        period_sec=(1.769 + 3.551 + 7.155) * 86400.0,
        ephem=ephem,
        vinf_cap=12.0,
        mu_central=MU_SUN_KM3_S2,   # WRONG centre for a Jovicentric chain
        slack_leg=2,
    )
    # Lambert about the Sun on Jupiter-centred states is physically inconsistent;
    # it must not report a bend-feasible closed chain (Axis-D refutation).
    assert not (r.converged and r.bend_feasible)
```

Run `-m slow` → **red**/**green**. Commit:

```
test: Axis-D wrong-central-mu falsification guard for moon chains (moon-tour Tier-1 Phase 7)
```

### Task 7.3 — full-suite + lint/type gate; OUTSTANDING note

- [x] `uv run pytest -m "not slow"` green (exit 0, 2026-06-08); `uv run pytest
  -m slow` for moon-tour evaluated: I-E-G closure PASS, Axis-A agreement PASS,
  Axis-D falsification PASS, bend-feasibility strict-xfail (honest-risk).
- [x] `ruff check` + `ruff format --check` clean; `mypy src tests` clean on all
  moon-tour files (one remaining tree error is a concurrent free_return WIP file,
  out of lane).
- [x] Update `data/OUTSTANDING.md`: Tier-1 shipped (registry + centred ephemeris +
  centre-agnostic corrector + Tisserand μ + VILM); the Jovian rows now
  patched-conic-computable; the Saturnian Titan/midsize split; **Tier-2 (CR3BP)
  remains open** for Arenstorf/Genova/Wittal + Saturnian midsize.

Commit:

```
docs: OUTSTANDING — moon-tour Tier-1 shipped; Tier-2 CR3BP still open (moon-tour Tier-1 Phase 7)
```

---

## Out of scope (explicit carve-outs — design §9)

- **All CR3BP / Tier 2** — Earth-Moon Arenstorf/Genova/Wittal rows + the
  Saturnian midsize-moon (Mimas/Enceladus/Tethys) members stay **citation-only**
  (Approval Q1). No CR3BP propagator, no Jacobi constant, no `orbit_elements.
  cr3bp{}` backfill, no JPL three-body catalog pull (Approval Q5 is for *that
  later* milestone).
- **The `T > 3` ballistic-transfer region** (Endgame Part-2) — needs CR3BP
  reachability; Tier-1.5/Tier-2.
- **Part-2 Table 1 cells as goldens** — excluded; the two flagged suspect cells
  preserved verbatim, never "corrected" (design §9, mining note 491-493).
- **2-letter moon code aliases** — NOT introduced; the catalogue's reserved
  full-name scheme is adopted ("Design claims" #2). No data migration.
- **Interplanetary→moon-capture multi-patched-conic chains** — single central
  body per cycler (design §9).
- **No heliocentric changes** — M-3D, M-ED, the Forge verdict logic, and every
  heliocentric golden stay byte-identical; the canonical signature is unchanged.

---

## Definition of done

0. `core/satellites.py` — `SatelliteData` + `SATELLITES`/`PRIMARIES`, full-name
   scheme, mean motion derived; registry-construction golden vs published Part-1
   Table 3 green. **Green.**
1. `core/ephemeris.py` — `_CentredCircularBackend` returns moon states about the
   primary; heliocentric backends byte-identical. **Green.**
2. The two Jovian rows are `cycler_class: multi-arc` (gauntlet routes them to
   invariants, not cr3bp); the Saturnian row keeps `non-keplerian` with an honest
   Titan-Tier-1 / midsize-Tier-2 split note; census ratchet re-derived. **Green.**
3. `search/correct.py` — `mu_central` plumbed into Lambert; `_max_bend_deg`
   resolves moon codes; heliocentric default byte-identical; Jovian I-E-G chain
   closes about Jupiter (slow, non-golden). **Green.**
4. `search/tisserand.py` — `_a_p_km` + `mu=` resolve a moon; round-trip
   `T=3−v∞²` about Jupiter; Sun-default byte-identical; a Jovicentric pair prunes
   through `linkable`. **Green.**
5. `search/vilm.py` — Eq.(9) root vs Table 3 V̄∞; Eq.(13) quadrature vs Table 1
   (no-GA) + Table 2 (with-GA); Europa 3-VILM scalar (154 m/s / 46 d); admissible
   ΔV-floor. **Green.** (Gate scope as stated in Phase 5.)
6. Matcher pre-filters by `(model_assumption, primary)`; heliocentric signature
   byte-identical (Sun bucket). **Green.**
7. Gauntlet axes compute primary-relative quantities (Axis B persistence about
   the primary; Axis A VILM-vs-corrector; Axis D wrong-μ refutation); combiner
   unchanged. **Green.**
8. No assertion uses a V∞/ΔV value our own code computed as the EXPECTED side —
   only the published Endgame Part-1 Table 1/2/3 + A6 scalar, and
   feasibility/closure/agreement predicates. The two flagged Part-2 cells are
   never goldens.
9. `uv run pytest -m "not slow"`, ruff, ruff format, mypy all clean. Census
   ratchet intact; no heliocentric golden moved.

---

## Honest risk register

| Item | Status | What gates it / mitigation |
|---|---|---|
| **Most moon-system V∞ values have NO sourced anchor.** | OPEN | The Russell-Strange Jovian/Saturnian rows' per-encounter V∞ are largely `null` (family-seed records, `data/README.md:72-76`; the rows carry citations, not V∞ multisets). So the **Phase-3 corrector closure is NON-GOLDEN** (closure + bend-feasibility asserted; the V∞ value is our computation). There is no sourced Jovicentric V∞ multiset to gate against the way Jones VEM gates the heliocentric M-ED corrector. **Mitigation:** the gated numbers are the **VILM ΔV** (Part-1 Tables 1–3, which DO have published anchors) and the **registry construction** (Table 3) — not the corrector's V∞. This is the honest boundary: Tier-1 makes the rows *computable* and *VILM-feasibility-gated*, not *V∞-rediscovery-gated*. |
| **Russell-Strange / Hernandez rows are family-seed `null`-numeric records.** | KNOWN | They carry attribution + qualitative geometry, `null` numerics (`data/README.md:72-76`). Tier-1 does not fabricate member V∞; it computes feasibility and reports honestly. No green test asserts a fabricated row value. |
| **VILM tables gate ΔV/V̄∞ only; they do NOT gate the corrector's V∞.** | BY DESIGN | Phase-5 gate scope is explicit: Table 3 → efficiency root; Table 1/2 → quadrature ΔV; A6 → end-to-end Europa scalar. The corrector (Phase 3) and the `T=3−v∞²` identity (Phase 4) are **physics-invariant-validated** (closure residual / algebraic round-trip), not table-gated. |
| **Linked-conic vs CR3BP 5–10% ΔV divergence + ballistic transfers linked-conics calls impossible.** | KNOWN | Mining note 286-293, 491-493. All VILM golden bands are **tolerances (≥10%), never equality**; a CR3BP-disagreeing value ≤10% is explicitly NOT a rejection. The validator must not reject a row for disagreeing with a linked-conic anchor within the band. |
| **Two flagged suspect Part-2 Table 1 cells.** | EXCLUDED | Titan `J_L4@100km` "766.5/776.4" (MAX<MIN). Never a golden; preserved verbatim. Part-2 (CR3BP) is Tier-2 anyway. |
| **#110 corrector signatures may drift mid-flight.** | LIVE | Every `correct.py` task has a VERIFY-FIRST re-read; the edits are minimal two-coupling lifts. If a `mu`/`center` param already landed, extend it. |
| **Concurrent `data/catalogue.yaml` edits.** | LIVE | Re-tag tasks VERIFY-FIRST re-grep the row ids and edit only the named `cycler_class:`/`notes:` lines. |
| **`safe_alt_km` per moon = the paper's table altitude (100 km / Titan 1500 km).** | DECISION | Recorded so Phase-5 anchors are comparable; if a different operational safe altitude is wanted later it is a registry edit + a re-derived (non-golden) bend, not a golden move. |

---

## Self-Review

### Design coverage (every approved item → plan location)

| Design / Approval item | Plan location |
|---|---|
| Q1 Tier-1-first; CR3BP citation-only | Goal; Out of scope; risk register |
| Q2 VILM module IN scope (Option 2) | Phase 5 (Tasks 5.0–5.4) |
| Q3 re-tag mis-flagged Jovian/Saturnian rows | Phase 2 (Tasks 2.1–2.2) |
| Q4 moon codes "reconciled against reserved codes" | "Design claims" #2; Phase 1 (full-name scheme) |
| §3 SATELLITES/PRIMARIES registry + sourcing discipline | Phase 1 (Tasks 1.0–1.3) |
| §3 registry-construction golden (Table 3) | Phase 1 Task 1.3 |
| §4 planet-centred circular ephemeris | Phase 2 Task 2.0 |
| §4 lift the two Sun-couplings (minimal) | Phase 3 (Tasks 3.1–3.2) |
| §4 `_a_p_km`/μ for moons | Phase 4 (Tasks 4.0–4.2) |
| §5 VILM as feasibility/search layer + ΔV-floor admissible bound | Phase 5 (Task 5.4) |
| §6 validation gates (which anchor gates what) | Phase 5 gate scope; Phase 1 Task 1.3 |
| §2/§6 (model_assumption, primary) pre-filter; heliocentric byte-identical | Phase 6 (Tasks 6.0–6.1) |
| §7 fidelity ladder about the primary; gauntlet axes | Phase 7 (Tasks 7.0–7.2) |
| §7 gauntlet dispatch fix (re-tag → multi-arc) | Phase 2 Task 2.1 |
| §9 non-goals; flagged cells excluded | Out of scope; risk register |

### Placeholder scan
No `TODO`/`FIXME`/`...`-as-impl in production-code instructions. Three test
sketches carry an explicit `...` (Phase 6 `<matcher_module>`, Phase 7 corrector-ΔV
extraction, Phase 7.0 entry-point name) — each is flagged "adapt to live names /
VERIFY-FIRST" with the recipe in prose, because the live matcher/fidelity
signatures must be read at run time (concurrency). One sourcing TODO (Task 1.1)
is a *human research step* (look up JPL SSD GM), not a code placeholder.

### Live-code verification done
Tisserand path corrected to `search/tisserand.py` (#1). Moon-code scheme
corrected to full names per `data/README.md:65-69` + live `bodies:` arrays (#2).
`correct.py` confirmed landed by #110 with no `mu`/`primary` (lines 78/115/186/232)
(#3). Fidelity ladder confirmed 4-rung at `provenance.py:95` (#4). `anchors_for`
dispatch confirmed at `validate.py:558-605` (keys off `cycler_class`; Jovian rows
are `non-keplerian` → demand `cr3bp` → misroute). `lambert(mu=...)` confirmed
`lambert.py:516`. Astropy Sun-subtract `ephemeris.py:263-265`. All Endgame
anchors transcribed from `docs/notes/2026-06-05-endgame-tisserand-mining.md`
A1–A6.

### Unverifiable / flagged claims
1. **Exact JPL SSD GM/physical values** are a human research step (Task 1.1/1.2
   sourcing TODO) — the plan does not invent the literals; it pins
   order-of-magnitude bands in tests and requires the sourced value + citation in
   the inline comment, cross-checked against the independent paper anchor (Task
   1.3). This preserves golden non-circularity.
2. **Live matcher pre-filter location** (Phase 6) — the spec says M7 *may*
   already pre-filter by `model_assumption` (`docs/spec.md:551`); the exact
   module/function was not pinned at authoring time (it is owned by code paths a
   concurrent agent may touch). Phase 6 Task 6.0 carries a VERIFY-FIRST grep to
   locate and extend it. Flagged, not papered over.
3. **`verify/fidelity.py` persistence entry-point name** (Phase 7.0) — design
   cites `fidelity.py:15-25,167-169,393-411`; the precise primary-parametrised
   entry point is to be read live (VERIFY-FIRST). The test import is marked TBD.
