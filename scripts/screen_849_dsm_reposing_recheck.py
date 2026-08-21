"""#849: re-run #388's canonical DSM lane under the CORRECTED g/G posing.

`dsm_descriptor_seed._descriptor_params` used to read a row's
``free_return_arcs`` POSITIONALLY (``g_tofs[0]``/``g_tofs[1]``) to decide which
arc is the non-designated (lowercase g) arc and which is the designated
(uppercase G, Mars-transit) arc -- exactly the defect class #820 found and
fixed for ``campaign_russell12.py::build_genome``. #849 fixed
``_descriptor_params`` to identify the designated arc by descriptor CASE
(:func:`cyclerfinder.search.descriptor.designated_arc_index`, the same
#794/#820-established semantics, now shared by both lanes) instead of array
position.

This script runs :func:`cyclerfinder.search.dsm_descriptor_seed.close_row_dsm`
-- "#388's canonical DSM determination path" per its own catalogue framing --
on the real DE440 ephemeris across every ``free_return_arcs``-bearing row,
under BOTH the pre-#849 (buggy, positional) posing and the corrected posing,
so the two can be compared directly. The pre-#849 posing is replicated
verbatim from the removed code (matches #830's own uncorrected
``dsm_388_recheck.json`` numbers byte-for-byte on the rows that file covers --
this script additionally covers the 5 descriptor-bearing rows #830 never ran).

NO catalogue writeback here (this is the compute half of the established
compute/adjudicate split -- see #820/#826, #822/#828, #827/#854, #839/#855).

Run: uv run python scripts/screen_849_dsm_reposing_recheck.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import cyclerfinder.search.dsm_descriptor_seed as dds
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog

_OUT = Path("data/found/849_dsm_reposing_recheck/dsm_reposing_recheck.json")


def _old_descriptor_params(
    row: dict[str, Any],
) -> tuple[float, float, float, float, float, tuple[str, ...]] | None:
    """Verbatim replica of the pre-#849 ``_descriptor_params`` (positional read).

    ``g_tofs[0]``/``g_tofs[1]`` = the first two ``free_return_arcs`` entries
    carrying a non-null ``tof_years``, in raw catalogue list order -- no check
    of descriptor case (designated vs generic) or ``arc_type`` (generic vs
    half-rev) at all.
    """
    aph = (row.get("orbit_elements") or {}).get("aphelion_au")
    vinf_list = {
        e["body"]: float(e["vinf_kms"])
        for e in (row.get("vinf_kms_at_encounters") or [])
        if e.get("body") is not None and e.get("vinf_kms") is not None
    }
    fra = row.get("free_return_arcs") or []
    g_tofs = [a.get("tof_years") for a in fra if a.get("tof_years") is not None]
    seq_str: str | None = row.get("sequence_canonical")
    if (
        len(g_tofs) < 2
        or aph is None
        or "E" not in vinf_list
        or "M" not in vinf_list
        or not seq_str
    ):
        return None
    sequence = tuple(seq_str.split("-"))
    if len(sequence) < 2:
        return None
    return (
        float(aph),
        float(g_tofs[0]),
        float(g_tofs[1]),
        float(vinf_list["E"]),
        float(vinf_list["M"]),
        sequence,
    )


def _run_one(row: dict[str, Any], ephem: Ephemeris, *, posing_fn: Any) -> dict[str, Any]:
    new_fn = dds._descriptor_params
    dds._descriptor_params = posing_fn
    try:
        seed = dds.seed_dsm_chain_from_descriptor(row)
        res = dds.close_row_dsm(row, ephem)
    finally:
        dds._descriptor_params = new_fn
    out: dict[str, Any] = {
        "seeded": seed is not None,
        "converged": res.converged,
        "anchor_match": res.anchor_match,
        "hyperbolic_impossible": res.hyperbolic_impossible,
        "max_residual_kms": (
            None if res.max_residual_kms != res.max_residual_kms else res.max_residual_kms
        ),
        "vinf_anchor_kms": (
            None if res.vinf_anchor_kms != res.vinf_anchor_kms else res.vinf_anchor_kms
        ),
        "vinf_per_encounter_kms": list(res.vinf_per_encounter_kms),
        "dv_dsm_kms": list(res.dv_dsm_kms),
    }
    if seed is not None:
        out["arc_a_au"] = seed.arc_a_au
        out["arc_e"] = seed.arc_e
        out["transit_branch"] = seed.transit_branch
    return out


def main() -> None:
    cat = load_catalog()
    rows = [e.raw for e in cat.entries if e.raw.get("free_return_arcs")]
    ephem = Ephemeris("astropy")  # real DE440
    t0 = time.time()
    records: list[dict[str, Any]] = []
    for row in rows:
        rid = row["id"]
        old = _run_one(row, ephem, posing_fn=_old_descriptor_params)
        new = _run_one(row, ephem, posing_fn=dds._descriptor_params)
        flipped = (old["seeded"], old["converged"], old["anchor_match"]) != (
            new["seeded"],
            new["converged"],
            new["anchor_match"],
        )
        rec = {
            "id": rid,
            "validation_level": row.get("validation_level", "V0"),
            "sourced_vinf": {
                e["body"]: e["vinf_kms"]
                for e in (row.get("vinf_kms_at_encounters") or [])
                if e.get("body") in ("E", "M")
            },
            "old_uncorrected_posing": old,
            "new_corrected_posing": new,
            "verdict_flip": flipped,
        }
        records.append(rec)

        def _fmt(r: dict[str, Any]) -> str:
            return (
                f"seeded={r['seeded']!s:5} conv={r['converged']!s:5} match={r['anchor_match']!s:5}"
            )

        print(
            f"[{time.time() - t0:6.0f}s] {rid:24s} OLD {_fmt(old)} | NEW {_fmt(new)} "
            f"{'*** FLIP ***' if flipped else ''}",
            flush=True,
        )

    _OUT.parent.mkdir(parents=True, exist_ok=True)
    _OUT.write_text(json.dumps(records, indent=2) + "\n")

    n_flips = sum(1 for r in records if r["verdict_flip"])
    print(f"\n=== summary ({time.time() - t0:.0f}s) ===")
    print(f"rows: {len(records)}  verdict flips: {n_flips}")
    print(f"wrote: {_OUT}")


if __name__ == "__main__":
    main()
