"""#607 -- triple/quadruple small-body multi-moon system symmetric-closure sweep.

`#605` shortlist item 2: small-body systems with 2-3 confirmed moons around
the SAME primary -- (87) Sylvia (Romulus, Remus), (130) Elektra (3 moons,
a quadruple system), (45) Eugenia (Petit-Prince + S/2004 (45) 1, both
CONFIRMED per JPL SBDB), and (216) Kleopatra (AlexHelios, CleoSelene) -- have
zero published cycler-analogue record (per `#605`'s literature scan) and have
never been swept by this project's own moon-tour machinery, which has so far
only ever been pointed at PLANET-centric systems (Uranus, Jupiter, Saturn,
Neptune). This script points the SAME machinery at small-body primaries
instead.

Method choice (per the #607 task spec: check whether #549's genome OR #563's
generalizes to a 3+-moon repeated sequence)
-----------------------------------------------------------------------------
`#549`'s real-binary `(k1,k2)` genome (`search/real_binary_kk_sweep.py`) is a
CR3BP construction: it treats the system as exactly TWO gravitating bodies
(primary + secondary) with the spacecraft as a massless third body orbiting
in that ROTATING two-body frame. This is fundamentally a two-primary model --
it has no notion of a third or fourth gravitating body at all, so it CANNOT
be generalized to a 3+-moon repeated sequence; there is no "add a moon" knob
in a restricted three-body problem. This was verified by reading
`real_binary_kk_sweep.py` in full (task #607 investigation), not assumed.

`#563`'s direct symmetric-closure construction (`enumerate_563_symmetric_
closures.py`, a patched-conic PRIMARY-centered moon-tour: spacecraft on
Keplerian Lambert arcs between circular-coplanar moons, gravity-assist bend
AT each moon) is exactly this project's existing multi-moon-tour method
(already applied to Uranus/Jupiter/Saturn/Neptune) -- a natural fit, since a
small-body primary + several much-less-massive moons is architecturally the
SAME hierarchy as a planet + its regular moons, just at a different mass/
length scale. `#600`'s 3-moon-chain generalization
(`enumerate_600_3moon_symmetric_closures.py`) extends this to the untested
Elektra case (3 moons, needs an anchor->flyby1->flyby2->anchor sequence).
Both scripts already accept an arbitrary `primary=`/moon-list, added by the
prior #571/#575 genericization passes -- this script is a thin driver that
supplies the SOURCED small-body registry entries (`core/satellites.py`,
added by this same task) and calls that existing, already-validated
construction verbatim. No new gate or construction logic here.

Sourcing discipline (per the #607 task spec)
-----------------------------------------------------------------------------
Every GM/sma value used here is sourced in `core/satellites.py`'s own
comments (with paper/author/year citations), not invented. Individual moon
masses are frequently UNMEASURED (too small to perturb their primary or each
other detectably) -- where this is the case, an ASSUMED-DENSITY mass is used
(flagged explicitly in the registry comments), matching the precedent
`#549` already set for Didymos-Dimorphos. See `core/satellites.py` for the
full per-body citation trail.

Lempo-Paha-Hiisi EXCLUDED (structural, not a data-sourcing failure)
-----------------------------------------------------------------------------
The #605 shortlist named the trans-Neptunian triple Lempo-Paha-Hiisi as a
candidate. Sourced research for this task (Benecchi, Noll, Grundy & Levison
2010, Icarus 207, 978; Correia 2018, Icarus 305, 250, arXiv:1710.08401) shows
Lempo and Hiisi (the INNER, close pair) have a mass ratio of only ~1.27:1 --
a genuine near-equal-mass binary, not a dominant-primary-plus-small-moon
hierarchy. Paha (much smaller) orbits the Lempo-Hiisi BARYCENTER, not a
fixed "primary." This violates the core physical assumption underlying the
`#563`/`#600` construction (a point-mass-dominant primary with negligible-
mass Keplerian moons about it) for exactly the reason `#549` treated real
near-equal-mass binaries (Patroclus-Menoetius, mu=0.438) via full CR3BP
rather than a patched-conic moon tour. Forcing this system through a
primary-centered patched-conic model would silently misrepresent the
dynamics, so it is excluded here on structural grounds, not added to
`core/satellites.py`, and not swept.

Irregular-primary shape caveat (per the #607 task spec item 3)
-----------------------------------------------------------------------------
All four included primaries are demonstrably non-spherical (Sylvia:
triaxial 374x248x194 km, shape-derived J2~0.024; Elektra: 262x205x164 km,
J2~0.16-0.18 -- SEVERE; Eugenia: oblate, ~186-215 km depending on method;
Kleopatra: famous bilobate "dog-bone," ~270x94x78 km). This script uses a
POINT-MASS primary gravity field throughout (the same idealization every
other primary in this registry uses) -- NOT a J2/multipole correction. This
is an explicit, flagged limitation (per the task's own allowed alternative:
"or an honest caveat that the point-mass model is being used anyway"), not a
silent assumption of validity. The moon-orbit-radius-to-primary-size ratios
are all >=1.8x the primary's own half-length (Kleopatra, the tightest case)
and up to 5-13x the primary's mean radius (Elektra) -- comparable to or
looser than the Uranus-Miranda ratio (~5.1x) this registry already treats as
point-mass, so this is not an unprecedented approximation, but it is
explicitly NOT validated against a full multipole gravity model (which the
source papers needed for precision orbit fits).

Discipline: NO catalogue writeback, NO V1-V4-strict gauntlet run here.
No `preflight_search` call -- like `#563`/`#600`, this is a bounded, direct
CONSTRUCTION (nothing searched for; every candidate is finite and
enumerated), not an open-ended budget-gated sweep; it is not a `run_*.py`
script and is not subject to that AST ratchet.

Run as::

    uv run python scripts/enumerate_607_smallbody_multimoon_symmetric_closures.py
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
# Also on sys.path so `scripts.X` (dotted) resolves standalone, matching how
# tests/scripts/*.py already address these modules -- needed because a couple
# of sibling scripts (this one included) are reachable BOTH ways depending on
# caller, and mypy needs a single, consistent module identity to avoid a
# "Source file found twice under different module names" error.
sys.path.insert(0, str(ROOT))

from scripts.enumerate_563_symmetric_closures import enumerate_direction  # noqa: E402
from scripts.enumerate_600_3moon_symmetric_closures import enumerate_sequence  # noqa: E402

DATA_DIR = ROOT / "data"
OUT_PATH = DATA_DIR / "enumerate_607_smallbody_multimoon_symmetric_closures.jsonl"

# (primary, moons, kind) -- kind "2moon" uses #563's anchor->flyby->anchor
# construction (both directions per pair); kind "3moon" uses #600's
# anchor->flyby1->flyby2->anchor chain (all ordered 3-permutations).
# Lempo-Paha-Hiisi deliberately absent -- see module docstring.
SYSTEMS: list[tuple[str, tuple[str, ...], str]] = [
    ("Sylvia", ("Romulus", "Remus"), "2moon"),
    ("Eugenia", ("PetitPrince", "EugeniaS2"), "2moon"),
    ("Kleopatra", ("AlexHelios", "CleoSelene"), "2moon"),
    ("Elektra", ("ElektraBeta", "ElektraGamma", "ElektraDelta"), "3moon"),
]


def run_2moon_system(primary: str, moons: tuple[str, ...]) -> list[dict[str, Any]]:
    directions: list[tuple[str, str]] = []
    for a, b in itertools.combinations(moons, 2):
        directions.append((a, b))
        directions.append((b, a))
    results = []
    for anchor, flyby in directions:
        res = enumerate_direction(anchor, flyby, primary=primary)
        results.append(res)
    return results


def run_3moon_system(primary: str, moons: tuple[str, ...]) -> list[dict[str, Any]]:
    sequences = list(itertools.permutations(moons, 3))
    results = []
    for anchor, flyby1, flyby2 in sequences:
        res = enumerate_sequence(anchor, flyby1, flyby2, primary=primary)
        results.append(res)
    return results


def main() -> int:
    t0 = time.time()
    print("[607] small-body multi-moon symmetric-closure sweep", flush=True)
    print(
        f"[607] systems: {[(p, m, k) for p, m, k in SYSTEMS]}",
        flush=True,
    )
    print(
        "[607] Lempo-Paha-Hiisi EXCLUDED (Lempo:Hiisi mass ratio ~1.27:1, a genuine "
        "near-equal-mass binary -- see module docstring)",
        flush=True,
    )

    all_out: list[dict[str, Any]] = []
    total_evaluated = 0
    total_passes = 0
    per_system_summary: list[dict[str, Any]] = []

    for primary, moons, kind in SYSTEMS:
        t_sys = time.time()
        if kind == "2moon":
            results = run_2moon_system(primary, moons)
        else:
            results = run_3moon_system(primary, moons)
        sys_evaluated = sum(r["n_evaluated"] for r in results)
        sys_passes = sum(r["n_all_gates_passed"] for r in results)
        total_evaluated += sys_evaluated
        total_passes += sys_passes
        elapsed_sys = time.time() - t_sys
        print(
            f"[607] {primary} ({kind}, moons={moons}): "
            f"{len(results)} directions/sequences, evaluated={sys_evaluated} "
            f"all_gates_pass={sys_passes} ({elapsed_sys:.1f}s)",
            flush=True,
        )
        for res in results:
            if kind == "2moon":
                label = f"{res['anchor']}-{res['flyby']}-{res['anchor']}"
            else:
                label = f"{res['anchor']}-{res['flyby1']}-{res['flyby2']}-{res['anchor']}"
            print(
                f"[607]   {label}: evaluated={res['n_evaluated']} "
                f"sub_gate={res['n_subgate_residual_only']} "
                f"all_gates_pass={res['n_all_gates_passed']}",
                flush=True,
            )
            all_out.append({"kind": kind, "primary": primary, **res})
        per_system_summary.append(
            {
                "primary": primary,
                "moons": list(moons),
                "kind": kind,
                "n_directions_or_sequences": len(results),
                "n_evaluated": sys_evaluated,
                "n_all_gates_passed": sys_passes,
                "elapsed_s": elapsed_sys,
            }
        )

    elapsed = time.time() - t0
    print(
        f"[607] DONE: {total_evaluated} candidates directly evaluated across "
        f"{len(SYSTEMS)} systems, {total_passes} pass ALL gates (residual+bend+DOP853) "
        f"({elapsed:.1f}s)",
        flush=True,
    )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#607 small-body multi-moon symmetric-closure sweep",
                    "systems": per_system_summary,
                    "excluded": [
                        {
                            "system": "Lempo-Paha-Hiisi",
                            "reason": (
                                "Lempo:Hiisi mass ratio ~1.27:1 (near-equal-mass binary, "
                                "not primary+test-particle-moon) -- structurally violates "
                                "the point-mass-primary construction; see module docstring"
                            ),
                        }
                    ],
                    "total_evaluated": total_evaluated,
                    "total_all_gates_passed": total_passes,
                    "elapsed_s": elapsed,
                }
            )
            + "\n"
        )
        for res in all_out:
            kind = res["kind"]
            if kind == "2moon":
                summary = {
                    "kind": "direction_summary",
                    "primary": res["primary"],
                    "anchor": res["anchor"],
                    "flyby": res["flyby"],
                    "t_syn_days": res["t_syn_days"],
                    "n_max": res["n_max"],
                    "n_evaluated": res["n_evaluated"],
                    "n_infeasible": res["n_infeasible"],
                    "n_subgate_residual_only": res["n_subgate_residual_only"],
                    "n_all_gates_passed": res["n_all_gates_passed"],
                }
            else:
                summary = {
                    "kind": "sequence_summary",
                    "primary": res["primary"],
                    "anchor": res["anchor"],
                    "flyby1": res["flyby1"],
                    "flyby2": res["flyby2"],
                    "sequence": res["sequence"],
                    "n_max_a": res["n_max_a"],
                    "n_max_b": res["n_max_b"],
                    "n_max_c": res["n_max_c"],
                    "n_evaluated": res["n_evaluated"],
                    "n_infeasible": res["n_infeasible"],
                    "n_subgate_residual_only": res["n_subgate_residual_only"],
                    "n_all_gates_passed": res["n_all_gates_passed"],
                }
            fh.write(json.dumps(summary) + "\n")
            for p in res["passes"]:
                fh.write(json.dumps({"kind": "pass", "primary": res["primary"], **p}) + "\n")
    print(f"[607] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
