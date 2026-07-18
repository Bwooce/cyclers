#!/usr/bin/env python3
"""#642: audit whether `#641`'s degenerate-equilibrium contamination finding
also affects `#608`'s/`#624`'s ORIGINAL Earth-Moon/cross-mu lift measurements.

`#641` (a Sun-Jupiter CR3BP periodic-orbit census) found that
``is_physically_sane`` had no check for degenerate Lagrange-point equilibria:
a fixed point of the CR3BP rotating-frame equations (velocity identically
zero) trivially satisfies ``correct_periodic``'s periodicity residual for
*any* period guess, so the corrector reports ``converged=True`` even though
nothing is actually orbiting. At Sun-Jupiter mu, this contaminated 609/614
(99%!) of "physically-sane converged" results.

This matters because `#624`'s cross-mu lift numbers (30x at mu=0.001, 3.5x at
Sun-Earth mu) are the specific evidence that got `#542` upgraded from
"answered" to "validated discovery lever worth productionizing". If those
numbers are inflated by silently counting fixed-point seeds as real converged
orbits, the actual lift -- and possibly the "lift genuinely transfers cross-mu"
conclusion itself -- could be weaker than believed.

**This script does NOT re-run `#608`'s or `#624`'s original (expensive)
searches.** Both tasks persisted their raw per-candidate refine results,
including the full converged ``state0`` (not just summary statistics):

  * ``data/found/608_generative_seed_poc/refine_results.jsonl``
  * ``data/found/624_cross_mu_transfer_pilot/refine_results.jsonl``

so the corrected numbers can be re-derived by re-filtering the EXACT SAME
saved converged states through the now-fixed
``cyclerfinder.ml.orbit_generative.is_physically_sane`` (which as of `#642`
rejects degenerate equilibria by default via
:func:`~cyclerfinder.ml.orbit_generative.is_degenerate_equilibrium`) --
deterministic, exact, and fast (no corrector calls at all).

**No catalogue writeback** -- this is a measurement-integrity audit of two
already-closed capability-evaluation tasks, not a discovery result.

**preflight_search() exemption**: this script re-reads two existing, already
-persisted result artifacts and re-filters them; there is no
``region_id``/``n_points`` sweep-region concept to preflight -- same category
as `#608`/`#614`/`#317`/`#624`/`#641` before it (see
``tests/scripts/test_scripts_call_preflight.py``'s ``_LEGACY_EXEMPT`` entry
for this file).
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cyclerfinder.ml.orbit_generative import (
    is_degenerate_equilibrium,
    is_physically_sane,
    lagrange_point_label,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent
_608_ARTIFACT = _REPO_ROOT / "data" / "found" / "608_generative_seed_poc" / "refine_results.jsonl"
_624_ARTIFACT = (
    _REPO_ROOT / "data" / "found" / "624_cross_mu_transfer_pilot" / "refine_results.jsonl"
)
_OUT_DIR = _REPO_ROOT / "data" / "found" / "642_equilibrium_contamination_audit"

# #608's Earth-Moon training mu; #624's two cross-mu targets.
_MU_EARTH_MOON = 0.01215058439469525
_MU_0_001 = 0.001
_MU_SUN_EARTH = 3.0034805950690393e-06


@dataclass(frozen=True)
class ArmResult:
    tag: str
    mu: float
    n: int
    n_converged: int
    n_sane_before: int  # PRE-#642 (equilibria counted as sane) -- original headline
    n_sane_after: int  # POST-#642 (equilibria rejected)
    n_equilibria_removed: int
    equilibrium_labels: dict[str, int]
    rate_before: float
    rate_after: float


def _load_records(path: Path, tags: set[str]) -> dict[str, list[dict[str, Any]]]:
    by_tag: dict[str, list[dict[str, Any]]] = {t: [] for t in tags}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("tag") in by_tag:
                by_tag[rec["tag"]].append(rec)
    return by_tag


def _analyze_arm(tag: str, mu: float, records: list[dict[str, Any]]) -> ArmResult:
    n = len(records)
    n_converged = 0
    n_sane_before = 0
    n_sane_after = 0
    equilibrium_labels: Counter[str] = Counter()
    for r in records:
        if not r.get("converged"):
            continue
        n_converged += 1
        state0, period, jacobi = r.get("state0"), r.get("period"), r.get("jacobi")
        if state0 is None or period is None or jacobi is None:
            continue
        sane_before = is_physically_sane(state0, period, jacobi, reject_degenerate_equilibria=False)
        if not sane_before:
            continue
        n_sane_before += 1
        if is_degenerate_equilibrium(state0):
            equilibrium_labels[lagrange_point_label(state0, mu)] += 1
        else:
            n_sane_after += 1
    return ArmResult(
        tag=tag,
        mu=mu,
        n=n,
        n_converged=n_converged,
        n_sane_before=n_sane_before,
        n_sane_after=n_sane_after,
        n_equilibria_removed=n_sane_before - n_sane_after,
        equilibrium_labels=dict(equilibrium_labels),
        rate_before=(n_sane_before / n if n else 0.0),
        rate_after=(n_sane_after / n if n else 0.0),
    )


def _ratio(numerator: float, denominator: float) -> float | None:
    """Lift ratio, honestly reported as ``None`` (not a fabricated number or a
    crash) when the denominator is zero -- an undefined ratio is a real,
    reportable outcome, not a bug to paper over.
    """
    if denominator <= 0.0:
        return None
    return numerator / denominator


def _report_pair(label: str, gen: ArmResult, base: ArmResult) -> dict[str, Any]:
    lift_before = _ratio(gen.rate_before, base.rate_before)
    lift_after = _ratio(gen.rate_after, base.rate_after)
    print(f"\n=== {label} ===")
    print(
        f"  generated: {gen.n_sane_before}/{gen.n} sane BEFORE -> "
        f"{gen.n_sane_after}/{gen.n} sane AFTER "
        f"({gen.n_equilibria_removed} removed as equilibria {gen.equilibrium_labels})"
    )
    print(
        f"  baseline:  {base.n_sane_before}/{base.n} sane BEFORE -> "
        f"{base.n_sane_after}/{base.n} sane AFTER "
        f"({base.n_equilibria_removed} removed as equilibria {base.equilibrium_labels})"
    )
    if lift_before is not None:
        print(f"  lift BEFORE (original headline): {lift_before:.2f}x")
    else:
        print("  lift BEFORE (original headline): undefined (0 baseline)")
    if lift_after is not None:
        print(f"  lift AFTER  (#642-corrected):    {lift_after:.2f}x")
    else:
        print("  lift AFTER  (#642-corrected):    undefined (0 baseline)")
    return {
        "label": label,
        "generated": asdict(gen),
        "baseline": asdict(base),
        "lift_before": lift_before,
        "lift_after": lift_after,
    }


def main() -> None:
    if not _608_ARTIFACT.exists() or not _624_ARTIFACT.exists():
        raise SystemExit(
            f"missing raw artifact(s): expected both {_608_ARTIFACT} and {_624_ARTIFACT} to "
            "exist (this audit re-derives numbers from #608's/#624's own saved converged "
            "states, never re-running their original searches)."
        )

    recs_608 = _load_records(_608_ARTIFACT, {"generated", "baseline_uniform"})
    recs_624 = _load_records(
        _624_ARTIFACT,
        {"mu=0.001_generated", "mu=0.001_baseline", "Sun-Earth_generated", "Sun-Earth_baseline"},
    )

    r608_gen = _analyze_arm("608_generated", _MU_EARTH_MOON, recs_608["generated"])
    r608_base = _analyze_arm("608_baseline", _MU_EARTH_MOON, recs_608["baseline_uniform"])
    r624_001_gen = _analyze_arm("624_mu0.001_generated", _MU_0_001, recs_624["mu=0.001_generated"])
    r624_001_base = _analyze_arm("624_mu0.001_baseline", _MU_0_001, recs_624["mu=0.001_baseline"])
    r624_se_gen = _analyze_arm(
        "624_sun_earth_generated", _MU_SUN_EARTH, recs_624["Sun-Earth_generated"]
    )
    r624_se_base = _analyze_arm(
        "624_sun_earth_baseline", _MU_SUN_EARTH, recs_624["Sun-Earth_baseline"]
    )

    print("#642: degenerate-equilibrium contamination audit of #608's/#624's raw saved states")
    print(f"  #608 artifact: {_608_ARTIFACT}")
    print(f"  #624 artifact: {_624_ARTIFACT}")

    results = [
        _report_pair("608 Earth-Moon in-distribution", r608_gen, r608_base),
        _report_pair("624 mu=0.001 cross-mu", r624_001_gen, r624_001_base),
        _report_pair("624 Sun-Earth cross-mu", r624_se_gen, r624_se_base),
    ]

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps({"task": "#642", "results": results}, indent=2))
    print(f"\nwrote {summary_path}")


if __name__ == "__main__":
    main()
