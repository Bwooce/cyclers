"""Re-verify past tulip JSONL rows under the #322 3D-topology gate.

Reads previously-emitted tulip-search JSONLs (from #281 / #283 / #313) and
re-evaluates each row's reported "tulip" claim using the FIXED gate from
:func:`cyclerfinder.genome.tulip.is_three_dimensional`. Emits one summary row
per re-checked candidate to ``data/tulip_topology_reverify_322.jsonl``.

This is a REGISTRY update, not a catalogue writeback. The catalogue
``data/catalogue.yaml`` is NOT modified.

Per the project memory ``feedback_bugfix_invalidates_past_searches``: every
correctness fix triggers a sweep + re-run of past searches computed under the
buggy code. Output schema:

    {
      "source_jsonl": "data/tulip_sweep_281.jsonl",
      "source_row_idx": 0,
      "system": "earth-moon",
      "np_target": 2,
      "original_converged": true,
      "original_reason": "ok",
      "original_z0": 0.1729,
      "max_abs_z": 0.183,
      "is_3d_under_fix": true,
      "verdict_under_fix": "3D tulip" | "planar Np-petal collapse" | ...,
      "downgraded_to_planar_collapse": false,
      "z_floor_nondim": 5e-3,
      "notes": "..."
    }
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.tulip import (
    TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
    is_three_dimensional,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"

# JSONLs to re-verify and the (primary, secondary) extraction helper.
SOURCES = [
    ("data/tulip_sweep_281.jsonl", "tulip_sweep_281"),
    ("data/tulip_higher_np_283.jsonl", "tulip_higher_np_283"),
    ("data/scan_313_mars_phobos.jsonl", "scan_313_mars_phobos"),
    ("data/scan_313_mars_deimos.jsonl", "scan_313_mars_deimos"),
    ("data/scan_313_sun_jupiter_io.jsonl", "scan_313_sun_jupiter_io"),
    ("data/scan_313_sun_jupiter_europa.jsonl", "scan_313_sun_jupiter_europa"),
]


_SYSTEM_TO_PAIR: dict[str, tuple[str, str]] = {
    "earth-moon": ("Earth", "Moon"),
    "mars-phobos": ("Mars", "Phobos"),
    "mars-deimos": ("Mars", "Deimos"),
    "jupiter-io": ("Jupiter", "Io"),
    "jupiter-europa": ("Jupiter", "Europa"),
    "jupiter-ganymede": ("Jupiter", "Ganymede"),
    "jupiter-callisto": ("Jupiter", "Callisto"),
    "saturn-enceladus": ("Saturn", "Enceladus"),
    "saturn-rhea": ("Saturn", "Rhea"),
    "saturn-titan": ("Saturn", "Titan"),
    "saturn-iapetus": ("Saturn", "Iapetus"),
    "uranus-titania": ("Uranus", "Titania"),
    "uranus-oberon": ("Uranus", "Oberon"),
    "neptune-triton": ("Neptune", "Triton"),
    "pluto-charon": ("Pluto", "Charon"),
    "sun-jupiter-io": ("Jupiter", "Io"),
    "sun-jupiter-europa": ("Jupiter", "Europa"),
}


def _row_primary_secondary(row: dict[str, Any]) -> tuple[str, str]:
    """Derive (primary, secondary) from a row, accepting both explicit fields
    and the ``system`` slug (e.g. ``"saturn-titan"`` -> ``("Saturn", "Titan")``).
    """
    if "primary" in row and "secondary" in row:
        return str(row["primary"]), str(row["secondary"])
    system_label = str(row.get("system", "")).lower()
    if system_label in _SYSTEM_TO_PAIR:
        return _SYSTEM_TO_PAIR[system_label]
    # Fallback: split on first hyphen and Title-case.
    parts = system_label.split("-")
    if len(parts) >= 2:
        return parts[0].title(), parts[-1].title()
    raise KeyError(f"could not infer primary/secondary from row: {row.get('system')}")


def _row_state0_period_system(row: dict[str, Any]) -> tuple[np.ndarray, float, cr3bp.CR3BPSystem]:
    """Reconstruct the (state0, period, system) triple from a JSONL row.

    All sources follow the same convention: perpendicular-crossing IC
    ``(x0, 0, z0, 0, ydot0, 0)``, period in nondim TU. The system pair is
    derived from explicit ``primary``/``secondary`` fields if present,
    otherwise from the ``system`` slug via :data:`_SYSTEM_TO_PAIR`.
    """
    primary, secondary = _row_primary_secondary(row)
    system = cr3bp.cr3bp_system(primary, secondary)
    x0 = float(row["x0"])
    z0 = float(row["z0"])
    ydot0 = float(row["ydot0"])
    period = float(row.get("T_nondim") or row.get("T_TU") or 0.0)
    state0 = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)
    return state0, period, system


_CLAIM_REASONS = {"ok", "direct_seed_match"}


def _row_is_a_claimed_tulip(row: dict[str, Any]) -> bool:
    """A row is "claimed tulip" iff:

      * ``converged == True`` (the search reported success), AND
      * the success reason is in :data:`_CLAIM_REASONS` (the success was
        attributed to a tulip topology match, NOT a fallback that the caller
        rejected anyway), AND
      * x0, z0, ydot0 are present (we need the IC to re-verify).

    Note: we DO NOT require ``n_obs == np_target`` here -- the #313
    scan_313 scripts mis-label rows (``np_target`` field carries the SEED
    family Np, not the petal target asked of find_tulip_at_system), so a
    row with ``n_obs=2, np_target=4, reason="direct_seed_match"`` IS a
    claimed tulip from the buggy gate's perspective.
    """
    if not bool(row.get("converged", False)):
        return False
    reason = str(row.get("reason") or row.get("notes") or "")
    # Match prefix only -- "ok" / "direct_seed_match" / sometimes prefixed
    # by other context.
    if not any(reason.startswith(r) or reason == r for r in _CLAIM_REASONS):
        return False
    return all(row.get(k) is not None for k in ("x0", "z0", "ydot0"))


def _system_label(row: dict[str, Any]) -> str:
    s = row.get("system")
    if isinstance(s, str):
        return s
    return f"{row.get('primary', '?')}-{row.get('secondary', '?')}".lower()


def reverify_jsonl(path: Path) -> list[dict[str, Any]]:
    """Re-verify all rows in a JSONL file; return a list of summary dicts."""
    if not path.exists():
        return []
    rows_out: list[dict[str, Any]] = []
    with path.open() as fh:
        for idx, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if not _row_is_a_claimed_tulip(row):
                # Not a claimed tulip; nothing to downgrade. Still log a stub
                # for full provenance.
                rows_out.append(
                    {
                        "source_jsonl": str(path.relative_to(REPO_ROOT)),
                        "source_row_idx": idx,
                        "system": _system_label(row),
                        "np_target": row.get("np_target"),
                        "original_converged": bool(row.get("converged", False)),
                        "original_reason": row.get("reason") or row.get("notes"),
                        "original_z0": row.get("z0"),
                        "max_abs_z": None,
                        "is_3d_under_fix": None,
                        "verdict_under_fix": "n/a (not a claimed tulip)",
                        "downgraded_to_planar_collapse": False,
                        "z_floor_nondim": TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
                        "notes": "skipped: row did not claim a tulip in the first place",
                    }
                )
                continue
            try:
                state0, period, system = _row_state0_period_system(row)
                is_3d, max_abs_z = is_three_dimensional(state0, period, system)
            except Exception as exc:
                rows_out.append(
                    {
                        "source_jsonl": str(path.relative_to(REPO_ROOT)),
                        "source_row_idx": idx,
                        "system": _system_label(row),
                        "np_target": row.get("np_target"),
                        "original_converged": True,
                        "original_reason": row.get("reason") or row.get("notes"),
                        "original_z0": row.get("z0"),
                        "max_abs_z": None,
                        "is_3d_under_fix": None,
                        "verdict_under_fix": "reverify_error",
                        "downgraded_to_planar_collapse": False,
                        "z_floor_nondim": TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
                        "notes": f"exception:{type(exc).__name__}:{exc}",
                    }
                )
                continue
            verdict = "3D tulip" if is_3d else "planar Np-petal collapse"
            rows_out.append(
                {
                    "source_jsonl": str(path.relative_to(REPO_ROOT)),
                    "source_row_idx": idx,
                    "system": _system_label(row),
                    "np_target": row.get("np_target"),
                    "original_converged": True,
                    "original_reason": row.get("reason") or row.get("notes"),
                    "original_z0": float(row["z0"]),
                    "max_abs_z": max_abs_z,
                    "is_3d_under_fix": is_3d,
                    "verdict_under_fix": verdict,
                    "downgraded_to_planar_collapse": (not is_3d),
                    "z_floor_nondim": TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
                    "notes": "",
                }
            )
    return rows_out


def main() -> int:
    out_rows: list[dict[str, Any]] = []
    for rel_path, label in SOURCES:
        path = REPO_ROOT / rel_path
        print(f"[reverify] {label}: reading {rel_path}")
        rows = reverify_jsonl(path)
        out_rows.extend(rows)
        n_claimed = sum(
            1
            for r in rows
            if r["is_3d_under_fix"] is not None
            or r["verdict_under_fix"] != "n/a (not a claimed tulip)"
        )
        n_downgraded = sum(1 for r in rows if r["downgraded_to_planar_collapse"])
        print(f"[reverify] {label}: claimed_tulip={n_claimed} downgraded={n_downgraded}")

    out_path = DATA_DIR / "tulip_topology_reverify_322.jsonl"
    with out_path.open("w") as fh:
        for r in out_rows:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
    print(f"[reverify] wrote {len(out_rows)} rows -> {out_path}")

    # Summary across all sources.
    print("\n=== Summary ===")
    by_source: dict[str, dict[str, int]] = {}
    for r in out_rows:
        s = r["source_jsonl"]
        if s not in by_source:
            by_source[s] = {"total": 0, "claimed": 0, "downgraded": 0, "still_3d": 0}
        by_source[s]["total"] += 1
        if r["is_3d_under_fix"] is not None:
            by_source[s]["claimed"] += 1
            if r["downgraded_to_planar_collapse"]:
                by_source[s]["downgraded"] += 1
            else:
                by_source[s]["still_3d"] += 1
    for src, counts in by_source.items():
        print(
            f"  {src}: total={counts['total']} "
            f"claimed_tulip={counts['claimed']} "
            f"still_3d={counts['still_3d']} "
            f"downgraded_to_planar={counts['downgraded']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
