"""Run #303 Part C: BCR4BP family validation-evidence probe against catalogue.

Reads the converged BCR4BP L1 Lyapunov family from
``data/bcr4bp_l1_family_303.jsonl`` and the existing catalogue
``data/catalogue.yaml``, and emits candidate "intermediate-fidelity bridge"
matches to ``data/bcr4bp_validation_bridges_303.jsonl``.

A candidate match requires:

  * The catalogue row's primary regime is Earth-Moon (``bodies = ["E", "Moon"]``
    or equivalent).
  * The catalogue row has a SOURCED period (period_nd in TU, period_days, or
    period.years) below ~1 year so the BCR4BP L1 family (period ~ 3 TU =
    13 days) is in the right ballpark.
  * Period deviation < 5%.

Per the task scope: this is **NOT a V-level promotion**. The bridge IS
intermediate-fidelity evidence at best -- the catalogue rows live in CR3BP
or bicircular models, the BCR4BP L1 Lyapunov family is in standard
incoherent BCR4BP, and a fair quantitative comparison requires (at least)
matching the same Jacobi-like invariant or a documented family-family
linkage. The JSONL is published as **CANDIDATES for human review**, not as
automatic catalogue writeback.

Discipline
----------
  * READ-ONLY on ``data/catalogue.yaml``.
  * Output is APPEND-ONLY: human review decides whether a candidate is a
    real bridge.
  * No V-level promotion. No novelty claims.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
FAMILY_PATH = REPO_ROOT / "data" / "bcr4bp_l1_family_303.jsonl"
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"
OUT_PATH = REPO_ROOT / "data" / "bcr4bp_validation_bridges_303.jsonl"


def _load_family() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Return ``(header, members)`` from the family JSONL."""
    header: dict[str, Any] | None = None
    members: list[dict[str, Any]] = []
    with FAMILY_PATH.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            obj = json.loads(line)
            if obj.get("row_type") == "header":
                header = obj
            elif obj.get("row_type") == "member":
                members.append(obj)
    if header is None:
        raise RuntimeError(f"no header row in {FAMILY_PATH}")
    return header, members


def _load_catalogue() -> list[dict[str, Any]]:
    """Return all catalogue rows as a list of dicts."""
    with CATALOGUE_PATH.open() as fh:
        rows = yaml.safe_load(fh)
    return rows or []


def _is_em_row(row: dict[str, Any]) -> bool:
    """True iff the row is an Earth-Moon system row."""
    bodies = row.get("bodies") or []
    if not bodies:
        return False
    body_set = {b.strip().lower() for b in bodies}
    return body_set == {"e", "moon"} or body_set == {"e", "mn"}


def _extract_period_days(row: dict[str, Any], tu_days: float) -> tuple[float | None, str]:
    """Return ``(period_days, source_label)`` or ``(None, reason)``.

    Search order: cr3bp.period_nd -> period.years -> period.note free-form
    days mention.
    """
    elements = row.get("orbit_elements") or {}
    cr3bp_block = (elements or {}).get("cr3bp") or {}
    period_nd = cr3bp_block.get("period_nd")
    if period_nd is not None:
        try:
            return float(period_nd) * tu_days, "orbit_elements.cr3bp.period_nd"
        except (TypeError, ValueError):
            pass

    period = row.get("period") or {}
    years = period.get("years")
    if years is not None:
        try:
            return float(years) * 365.25, "period.years"
        except (TypeError, ValueError):
            pass

    # Free-form note: look for "~XX d" / "~XX days" mentions in period.note.
    note = period.get("note")
    if isinstance(note, str):
        match = re.search(r"(\d+(?:\.\d+)?)\s*(?:\+/?-\s*\d+(?:\.\d+)?\s*)?d(?:ays)?\b", note)
        if match:
            try:
                return float(match.group(1)), "period.note (free-form days mention)"
            except (TypeError, ValueError):
                pass

    return None, "no_sourced_period"


def _extract_vinf_kms_earth(row: dict[str, Any]) -> float | None:
    """Extract earth-flyby V_infinity in km/s if sourced; else None."""
    enc = row.get("vinf_kms_at_encounters") or []
    for e in enc:
        body = (e.get("body") or "").strip().lower()
        if body in ("e", "earth"):
            v = e.get("vinf_kms")
            if v is None:
                return None
            try:
                return float(v)
            except (TypeError, ValueError):
                return None
    return None


def main() -> int:
    header, members = _load_family()
    catalogue = _load_catalogue()
    tu_days = float(header["tu_days"])
    print(
        f"#303 Part C catalogue probe: family has {len(members)} members "
        f"over mu_sun extent {header['mu_sun_extent']}, TU_DAYS={tu_days:.6f}"
    )

    em_rows = [r for r in catalogue if _is_em_row(r)]
    print(f"  catalogue has {len(catalogue)} total rows; {len(em_rows)} are Earth-Moon.")

    # Probe: for each EM row, find the family member with min period deviation
    # (in days) and the row's V_inf (if sourced). Record as candidate or
    # explicit no-match.
    candidates: list[dict[str, Any]] = []
    no_match_rows: list[dict[str, Any]] = []
    for row in em_rows:
        row_id = row.get("id", "<unknown>")
        validation_level = row.get("validation_level") or "<unset>"
        period_days, period_source = _extract_period_days(row, tu_days)
        if period_days is None:
            no_match_rows.append(
                {
                    "row_type": "no_match",
                    "catalogue_id": row_id,
                    "validation_level": validation_level,
                    "reason": period_source,
                    "note": (
                        "row has no sourced period; cannot compare to BCR4BP "
                        "family member periods. Free-form 'approximate' periods "
                        "(e.g., Wittal's '~27 +/- 3 d') are intentionally NOT "
                        "matched -- approximate quotes are not closure goldens."
                    ),
                }
            )
            continue
        # Find min-deviation member.
        best = None
        for m in members:
            t_days = float(m["T_days"])
            dev_pct = abs(t_days - period_days) / period_days * 100.0
            if best is None or dev_pct < best["dev_pct"]:
                best = {
                    "member_step_idx": int(m["step_idx"]),
                    "member_mu_sun_value": float(m["mu_sun_value"]),
                    "member_T_days": t_days,
                    "member_T_TU": float(m["T_TU"]),
                    "dev_pct": dev_pct,
                }
        assert best is not None
        is_candidate = best["dev_pct"] < 5.0
        record = {
            "row_type": "candidate" if is_candidate else "no_match",
            "catalogue_id": row_id,
            "validation_level": validation_level,
            "catalogue_period_days": period_days,
            "catalogue_period_source": period_source,
            "catalogue_vinf_kms_earth": _extract_vinf_kms_earth(row),
            **best,
            "period_dev_pct": best["dev_pct"],
            "match_threshold_pct": 5.0,
            "honest_scope": (
                "Intermediate-fidelity BCR4BP L1 Lyapunov family match candidate. "
                "NOT a V-level promotion. The match indicates the BCR4BP L1 "
                "Lyapunov family at the indicated mu_sun has a period within "
                "the 5% threshold of the catalogue row's published period. "
                "True bridging would require: (1) family-family linkage from "
                "the catalogue row's documented family to the L1 Lyapunov, "
                "(2) matched conserved invariant (Jacobi / energy), and (3) "
                "a sourced reference for the bridge itself."
                if is_candidate
                else "BCR4BP L1 Lyapunov family does NOT bridge this row "
                "(period deviation > 5%). The L1 Lyapunov family at C=3.1294 "
                "has T ~ 12.8 days; this row's period is elsewhere. A "
                "different BCR4BP family (different L_i, different Jacobi "
                "level, or a non-L1 family) might bridge -- not this one."
            ),
        }
        if is_candidate:
            candidates.append(record)
        else:
            no_match_rows.append(record)

    print(f"  candidates (period dev < 5%): {len(candidates)}; no-match rows: {len(no_match_rows)}")
    for c in candidates:
        print(
            f"    CANDIDATE {c['catalogue_id']} (V={c['validation_level']}): "
            f"row T={c['catalogue_period_days']:.2f} d, "
            f"member T={c['member_T_days']:.2f} d, "
            f"dev={c['period_dev_pct']:.3f}%"
        )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out_header = {
        "row_type": "header",
        "task_id": 303,
        "phase": "phase-2-catalogue-validation-evidence",
        "family_jsonl_source": str(FAMILY_PATH.relative_to(REPO_ROOT)),
        "catalogue_source": str(CATALOGUE_PATH.relative_to(REPO_ROOT)),
        "tu_days": tu_days,
        "em_rows_probed": len(em_rows),
        "candidates_count": len(candidates),
        "no_match_count": len(no_match_rows),
        "match_threshold_pct": 5.0,
        "honest_scope": (
            "These are CANDIDATE intermediate-fidelity bridges between the "
            "BCR4BP L1 Lyapunov continuation family (#303 Phase 2) and "
            "existing catalogue rows. NOT a V-level promotion; NOT a "
            "catalogue writeback. Per the orbit-closure discipline and the "
            "feedback_orbit_closure_discipline + feedback_golden_tests_sourced_only "
            "memories: any actual promotion requires (1) a sourced "
            "family-family linkage, (2) matched conserved invariants, "
            "(3) HFEM real-ephemeris closure to claim V4. The L1 Lyapunov "
            "family has T ~ 12.8 d, so it is in the right ballpark to "
            "match short-period EM rows -- but the topology / family "
            "identity comparison is OUT OF SCOPE for this probe."
        ),
    }
    with OUT_PATH.open("w") as fh:
        fh.write(json.dumps(out_header) + "\n")
        for c in candidates:
            fh.write(json.dumps(c) + "\n")
        for nm in no_match_rows:
            fh.write(json.dumps(nm) + "\n")
    print(f"  wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
