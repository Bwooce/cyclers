"""M7 Phase 2 coverage scan (#423) — which catalogue cycler rows can M7 measure?

Runs verify_real_closure(compute_tcm=True, n_cycles=2) over every cycler-type
catalogue row that carries a priority_date, recording whether the real-eph cycler
CONSTRUCTS (full leg data), CLOSES (M6b), and what horizon TCM M7 returns
(finite / inf-diverged). Prints a per-row table + summary. Read-only; no writeback.
"""

from __future__ import annotations

import warnings
from datetime import UTC, datetime

import yaml

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import CATALOGUE_PATH
from cyclerfinder.verify.real_closure import verify_real_closure

_ALDRIN_DEFAULT = datetime(2002, 8, 14, tzinfo=UTC)


def main() -> None:
    warnings.simplefilter("ignore")
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    if isinstance(rows, dict):
        rows = rows.get("cyclers") or rows.get("entries") or list(rows.values())
    ephem = Ephemeris("astropy")

    cyc_rows = [r for r in rows if isinstance(r, dict) and (r.get("type") or "cycler") == "cycler"]
    print(f"{datetime.now(UTC).isoformat()}  scanning {len(cyc_rows)} cycler-type rows\n")
    print(f"{'id':40s} {'vlevel':6s} {'v3_status':22s} {'horizon_mps':>12s}")
    print("-" * 84)

    n_measured = n_diverged = n_skip = 0
    for r in cyc_rows:
        cid = str(r.get("id", "?"))
        raw = r.get("raw") or {}
        vlevel = str(r.get("validation_level") or raw.get("validation_level") or "")
        pri = r.get("priority_date")
        pri_dt = None
        if isinstance(pri, str):
            try:
                pri_dt = datetime.fromisoformat(pri.replace("Z", "+00:00"))
                if pri_dt.tzinfo is None:
                    pri_dt = pri_dt.replace(tzinfo=UTC)
            except ValueError:
                pri_dt = None
        try:
            res = verify_real_closure(
                r,
                2,
                ephem,
                cycler_id=cid,
                signature_priority_date=pri_dt or _ALDRIN_DEFAULT,
                compute_tcm=True,
            )
        except Exception as e:
            print(f"{cid:40.40s} {vlevel:6s} {'ERROR':22s} {str(e)[:30]}")
            n_skip += 1
            continue
        h = res.horizon_tcm_mps
        hs = "inf" if h == float("inf") else f"{h:.1f}"
        print(f"{cid:40.40s} {vlevel:6s} {res.v3_status:22s} {hs:>12s}")
        if res.v3_status.startswith("v3-construction") or res.v3_status == "v3-no-real-window":
            n_skip += 1
        elif h == float("inf"):
            n_diverged += 1
        else:
            n_measured += 1

    print("-" * 84)
    print(
        f"{datetime.now(UTC).isoformat()}  DONE: "
        f"measured(finite TCM)={n_measured}  diverged(inf)={n_diverged}  "
        f"skip(construct/window/error)={n_skip}"
    )


if __name__ == "__main__":
    main()
