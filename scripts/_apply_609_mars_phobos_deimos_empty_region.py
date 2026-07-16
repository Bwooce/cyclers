"""#609 -- Mars Phobos-Deimos symmetric-closure empty-region stamp.

One-shot writer (prefixed `_` so it is treated as scratch and left untouched
by the harness, matching `_apply_571_titan_smallmoon_empty_regions.py`'s own
convention): appends ONE EmptyRegionReport to `data/empty_regions.jsonl` for
the Mars Phobos-Deimos pair, per `#605`'s shortlist item 4 ("hierarchical
cycler-of-cyclers", piloted via Mars Phobos-Deimos as the cheapest possible
target moon-cycler -- never swept before this task, a genuine roster gap).

This is NOT an analytic pre-check (contrast `_apply_571_...py`, which skipped
the sweep entirely): the real, already-genericized #563 direct
symmetric-closure enumeration
(`scripts/enumerate_563_symmetric_closures.py --primary Mars --moons
Phobos,Deimos`) was actually RUN, producing
`data/enumerate_609_mars_phobos_deimos_symmetric_closures.jsonl`. This script
only re-derives the summary statistics from that run (+ a live re-check of
the sub-gate survivors' own flyby-encounter bend, since the #563 script does
not persist per-candidate bend for residual-only sub-gate misses) and writes
the registry stamp -- no new physics, no new sweep.

Result: 0/512 candidates pass ALL gates (residual + #324 physical max-bend +
DOP853 cross-check) in EITHER direction (Phobos-Deimos-Phobos,
Deimos-Phobos-Deimos). 52/512 clear the residual gate alone (both bodies are
so close together and near-commensurate that the kinematic Lambert closure is
essentially exact, residual ~1e-13-1e-15 km/s) but EVERY one fails the #324
physical bend gate: the flyby moon's own achievable bend tops out at 0.0159
deg across all 52 sub-gate survivors, over 300x short of the
DEFAULT_MIN_USEFUL_BEND_DEG=5.0 deg floor. Same failure mode as #571
(Saturn-Titan-{Mimas,Enceladus,Tethys,Dione}) and #599
(Neptune-Triton-Proteus): a moon too low-GM to bend a flyby usefully -- here
BOTH bodies in the pair are undersized (Phobos GM 7.09e-4, Deimos GM
9.62e-5 km^3/s^2; core/satellites.py's own sourcing note already flags both
as "negligible bending").

Consequence for #609's hierarchical-cycler-of-cyclers pilot: with no
gate-passing Phobos-Deimos moon-cycler to compose with, the phase-matching
step (checking commensurability against a catalogued Earth-Mars heliocentric
cycler) cannot proceed -- there is no moon-cycler period to phase-match
against. This is reported as the #609 deliverable itself (a legitimate,
bounded negative), not a partial result.

Run as::

    uv run python scripts/_apply_609_mars_phobos_deimos_empty_region.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.data.empty_regions import (  # noqa: E402
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.search.physical_sanity import DEFAULT_MIN_USEFUL_BEND_DEG  # noqa: E402

SWEEP_DATA_PATH = ROOT / "data" / "enumerate_609_mars_phobos_deimos_symmetric_closures.jsonl"

git_sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()

MU_MARS = PRIMARIES["Mars"]
PHOBOS = SATELLITES["Phobos"]
DEIMOS = SATELLITES["Deimos"]

# Loaded from the actual #563-generic run (scripts/enumerate_563_symmetric_closures.py
# --primary Mars --moons Phobos,Deimos), not re-derived here.
with SWEEP_DATA_PATH.open(encoding="utf-8") as fh:
    lines = [json.loads(line) for line in fh]
meta = lines[0]
assert meta.get("_meta")
direction_summaries = [rec for rec in lines if rec.get("kind") == "direction_summary"]
assert len(direction_summaries) == 2, direction_summaries

# Per-candidate flyby-own bend across all 52 residual-sub-gate survivors
# (recomputed live from the same #558/#563 gate machinery the real run used;
# see module docstring -- the #563 JSONL only persists ALL-GATE passes, so a
# fresh gate_candidate() call over the same finite grid is needed to report
# the sub-gate survivors' own bend numbers, not a new sweep).
sys.path.insert(0, str(ROOT / "scripts"))
import itertools  # noqa: E402
import math  # noqa: E402

from refine_562_commensurability import synodic_period_days  # noqa: E402
from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    GATE_RESIDUAL_KMS,
    gate_candidate,
    residual_at_point,
)

from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day  # noqa: E402

N_REV_VALUES = range(4)
REL_OFFSETS_DEG = (0.0, 180.0)

max_bend_seen = 0.0
n_subgate_total = 0
for anchor, flyby in (("Phobos", "Deimos"), ("Deimos", "Phobos")):
    sat_a = SATELLITES[anchor]
    sat_b = SATELLITES[flyby]
    t_syn = synodic_period_days(MU_MARS, sat_a.sma_km, sat_b.sma_km, opposite_sense=False)
    n_a = _mean_motion_rad_day(MU_MARS, sat_a.sma_km)
    n_b = _mean_motion_rad_day(MU_MARS, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    sqrt_papb = math.sqrt(p_a * p_b)
    tof_max = 3.0 * sqrt_papb
    n_max = math.floor(2.0 * tof_max / t_syn)
    for n in range(1, n_max + 1):
        tof_days = n * t_syn / 2.0
        tof_scale = tof_days / sqrt_papb
        for n0, n1 in itertools.product(N_REV_VALUES, N_REV_VALUES):
            for rel in REL_OFFSETS_DEG:
                pt = residual_at_point(
                    anchor,
                    flyby,
                    rel_offset_deg=rel,
                    tof_scale=tof_scale,
                    n_rev=(n0, n1),
                    primary="Mars",
                )
                if pt is None or pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                    continue
                n_subgate_total += 1
                gated = gate_candidate(anchor, flyby, pt, primary="Mars")
                max_bend_seen = max(max_bend_seen, gated["max_bend_deg_per_encounter"][1])

assert n_subgate_total == sum(d["n_subgate_residual_only"] for d in direction_summaries)
assert meta["total_all_gates_passed"] == 0

t_syn_days = direction_summaries[0]["t_syn_days"]

method = MethodCapability(
    genome="#563 direct symmetric-closure enumeration (rel_offset in {0,180 deg}, "
    "exactly-commensurate tof = n*T_syn/2), genericized (#575/#599) to a primary= "
    "parameter -- pointed at Mars/Phobos/Deimos for the first time (#609; never "
    "swept in this project before, a genuine roster gap per #605).",
    corrector="#558/#563 residual_at_point Lambert-closure construction, gated by the "
    "two-sided #324 physical max-bend gate (search/physical_sanity.py) and a DOP853 "
    "cross-check",
    capability_tags=frozenset(
        {
            "ballistic",
            "circular",
            "coplanar",
            "patched-conic",
            "symmetric-closure-construction",
        }
    ),
    git_sha=git_sha,
)

report = EmptyRegionReport(
    region_id="mars-phobos-deimos-symmetric-closure-609-2026-07-16",
    family="Mars Phobos-Deimos 2-moon (anchor-flyby-anchor) symmetric-closure "
    "extension of the #563 direct-construction method -- first application at "
    "Mars, first pair where BOTH moons are individually undersized",
    centre="Mars",
    topologies=(
        {
            "sequence": ["Phobos", "Deimos", "Phobos"],
            "note": "both directions (Phobos-anchor and Deimos-anchor) enumerated",
        },
    ),
    method_capability=method,
    search_extent={
        "n_directions_total": 2,
        "n_rev_combos": 16,
        "rel_offsets_deg": list(REL_OFFSETS_DEG),
        "tof_scale_max_bound": 3.0,
        "points_total": meta["total_evaluated"],
        "ephem_model": "circular-coplanar idealized Mars model (#563 construction)",
        "center": "Mars",
        "t_syn_days": t_syn_days,
        "gate_residual_kms": GATE_RESIDUAL_KMS,
        "n_subgate_residual_only_total": n_subgate_total,
        "phobos_gm_km3_s2": PHOBOS.mu_km3_s2,
        "deimos_gm_km3_s2": DEIMOS.mu_km3_s2,
    },
    prune_gates=(
        "#563 symmetric-closure construction (residual_at_point Lambert closure)",
        f"residual gate ({GATE_RESIDUAL_KMS} km/s) -- {n_subgate_total}/"
        f"{meta['total_evaluated']} candidates passed this gate",
        "two-sided #324 physical max-bend gate "
        "(search/physical_sanity.py::candidate_passes_physical_gate) -- ALL "
        f"{n_subgate_total} sub-gate survivors fail here",
    ),
    result={
        "n_evaluated_total": meta["total_evaluated"],
        "n_subgate_residual_only_passed": n_subgate_total,
        "n_all_gates_passed": meta["total_all_gates_passed"],
        "max_flyby_own_bend_deg_across_subgate_survivors": round(max_bend_seen, 4),
        "min_useful_bend_deg_required": DEFAULT_MIN_USEFUL_BEND_DEG,
    },
    verdict=f"EMPTY (method-conditional, physically-grounded): 0/{meta['total_evaluated']} "
    "candidates pass all gates. Both Phobos (GM 7.087e-4 km^3/s^2) and Deimos (GM "
    "9.62e-5 km^3/s^2) are far too low-mass to deliver a useful gravity-assist bend "
    f"-- max achieved bend across all {n_subgate_total} residual-gate survivors is only "
    f"{round(max_bend_seen, 4)} deg, ~300x under the {DEFAULT_MIN_USEFUL_BEND_DEG} deg "
    "usefulness floor.",
    interpretation=(
        "Same failure mode that excluded Miranda from the Uranian #563 census and "
        "Proteus from the #599 Neptune Triton-Proteus census (both single-sided: one "
        "undersized moon out of a pair). The Mars Phobos-Deimos pair is a DOUBLY "
        "undersized case -- BOTH bodies in the pair are far below any useful-bend "
        "threshold (core/satellites.py's own sourcing note already flags Phobos/Deimos "
        "GM < 1 km^3/s^2 as making them 'poor gravity-assist / cycler bodies', "
        "confirmed here quantitatively rather than just asserted). The kinematic "
        "Lambert closure itself is essentially exact at every commensurate tof (residual "
        "~1e-13 to 1e-15 km/s, since the two moons' orbits are close together and the "
        "synodic period is short, ~0.427 days) -- geometrically the loop closes almost "
        "trivially -- but this is irrelevant to cycler-usefulness without a real bend at "
        "the flyby node. No #609 hierarchical phase-matching pilot can proceed against "
        "this pair: there is no gate-passing Phobos-Deimos moon-cycler period to "
        "phase-match a heliocentric Earth-Mars cycler's Mars-arrival epoch against. Per "
        "[[project_negative_results_registry]], this negative is conditional on this "
        "project's own two-sided #324 gate policy (every node must contribute useful "
        "bending) and on the circular-coplanar patched-conic idealization -- a genuinely "
        "different architecture (e.g. Phobos/Deimos as PASSIVE science targets riding a "
        "third, more massive perturbation source, or a full 3-body low-thrust/electric "
        "design not gated on ballistic bend at all) is a structurally different, "
        "uncovered case, not a contradiction of this finding."
    ),
    source_anchors="none (this is a search-methodology negative, not a literature-corpus "
    "claim; literature_check.py's own Wallace-Phobos CR3BP-rendezvous anchor already notes "
    "no published Sun-Mars-Phobos BCR4BP or Phobos-Deimos repeated-flyby cycler exists)",
    run={
        "date": "2026-07-16",
        "task": 609,
        "data_files": [str(SWEEP_DATA_PATH.relative_to(ROOT))],
        "note": "#605 shortlist item 4 pilot (hierarchical cycler-of-cyclers via Mars "
        "Phobos-Deimos): step 1 of the pilot (does a gate-passing Phobos-Deimos "
        "symmetric-closure moon-cycler exist at all) returns a clean negative, so the "
        "hierarchical phase-matching step (step 2) cannot be attempted. No catalogue.yaml "
        "edit.",
    },
)

append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
print(f"[609] appended region_id={report.region_id!r} to {DEFAULT_EMPTY_REGIONS_PATH}")
