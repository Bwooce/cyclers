"""#338 Part B -- boundary analysis on the annual V4-strict sweep.

Reads ``data/silver_327_v4_strict_annual_sweep_338.jsonl`` (Part A) and
asks the headline question:

    Is the PASS pattern across the 100 launch epochs CYCLIC (suggesting
    a moon-phase commensurability the SILVER is locked to), or
    IRREGULAR / CHAOTIC (suggesting the bounded-drift signature is
    fragile and not mission-useful)?

What we compute
---------------
1. PASS/FAIL totals + per-decade breakdown.
2. FAIL-run gap statistics + clustering (are FAILs clumped or evenly
   spread?).
3. Autocorrelation of ``pass_v4_strict`` at lags 1..50 (look for any
   integer-year periodicity in the binary PASS series).
4. Autocorrelation of the continuous ``drift_agreement_kms_vs_v3``
   signal (finer than the binary PASS series).
5. Best-period candidates against physical moon-phase commensurabilities:

   * Umbriel orbital period ~4.144 days
   * Oberon orbital period ~13.46 days
   * Umbriel-Oberon synodic period ~5.99 days
   * Uranus heliocentric orbital period ~84.02 yr

   The first three are SUB-YEAR and the annual sweep cannot resolve
   them; we document the aliasing. The 84 yr period IS resolvable but
   our window is 100 yr -- only ~1.2 cycles, so the autocorrelation has
   limited statistical power.

6. Boundary verdict:

   * CYCLIC if a clear period emerges in either autocorrelation
     (peak >= 0.5 at a non-zero lag and lag is not trivially close
     to the window edge).
   * IRREGULAR if no clear period AND FAILs are clustered near the
     kernel-coverage edge (suggesting kernel-edge effects, not
     intrinsic resonance failure).
   * EFFECTIVELY_CYCLIC (special case) if the FAIL fraction is small
     AND fail clustering is confined to a known boundary (e.g. the
     URA111 kernel's 2099 endpoint).

Output
------
``data/silver_327_v4_strict_boundary_338.jsonl`` -- one row per analysis
section (PASS/FAIL totals, per-decade, autocorrelation, verdict).

NO catalogue writeback.

Run as::

    uv run python scripts/analyze_338_boundary.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

IN_JSONL = ROOT / "data" / "silver_327_v4_strict_annual_sweep_338.jsonl"
OUT_JSONL = ROOT / "data" / "silver_327_v4_strict_boundary_338.jsonl"

# Periodicity peak threshold for the CYCLIC verdict. Pearson autocorrelation
# at a non-trivial lag must exceed this for us to call the pattern cyclic.
AUTOCORR_PEAK_THRESHOLD = 0.50

# Don't trust autocorrelation peaks at lags >= LAG_FRACTION_MAX * N
# (statistical power drops as the lag approaches the window size).
LAG_FRACTION_MAX = 0.5

# Kernel-edge guard: FAILs in the last KERNEL_EDGE_YEARS of the sweep
# (where kernel coverage runs out at 2099) are flagged as kernel-edge
# artifacts and downweighted in the verdict.
KERNEL_EDGE_YEAR = 2099
KERNEL_EDGE_YEARS = 15  # i.e. 2085+ is "near the edge"

# Physical moon-phase commensurability candidates (days; reference values
# from Murray-Dermott Table A.7 + standard derivations).
UMBRIEL_PERIOD_DAYS = 4.144
OBERON_PERIOD_DAYS = 13.463
UMBRIEL_OBERON_SYNODIC_DAYS = 5.987  # = 1 / (1/T_U - 1/T_O)
URANUS_ORBITAL_YEARS = 84.02

# Acceptable launch-epoch sub-band: must have >= MIN_PASS_FRACTION PASSes
# AND the recommended launch_epoch must be in the middle of a stretch
# of PASSes (defensive scheduling).
MIN_PASS_FRACTION = 0.85


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _read_sweep() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Read Part A's JSONL -> (per-year rows, header meta)."""
    rows: list[dict[str, Any]] = []
    header: dict[str, Any] = {}
    with IN_JSONL.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("_meta") and "task" in r and not header:
                header = r
            elif r.get("kind") == "annual_sweep_row":
                rows.append(r)
    rows.sort(key=lambda r: r["year"])
    return rows, header


def _decade_breakdown(rows: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    """Group PASS/FAIL by decade."""
    out: dict[str, dict[str, int]] = {}
    for r in rows:
        dec = f"{(r['year'] // 10) * 10}s"
        if dec not in out:
            out[dec] = {"pass": 0, "fail": 0}
        if r["passes_v4_strict"]:
            out[dec]["pass"] += 1
        else:
            out[dec]["fail"] += 1
    return out


def _fail_clustering(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Are FAILs clustered? Compute fail-run gap stats."""
    fails = [r["year"] for r in rows if not r["passes_v4_strict"]]
    if len(fails) < 2:
        return {
            "n_fails": len(fails),
            "fail_years": fails,
            "gaps_yr": [],
            "gap_mean_yr": None,
            "gap_min_yr": None,
            "gap_max_yr": None,
            "fail_year_min": fails[0] if fails else None,
            "fail_year_max": fails[-1] if fails else None,
        }
    gaps = [fails[i + 1] - fails[i] for i in range(len(fails) - 1)]
    return {
        "n_fails": len(fails),
        "fail_years": fails,
        "gaps_yr": gaps,
        "gap_mean_yr": float(np.mean(gaps)),
        "gap_min_yr": int(min(gaps)),
        "gap_max_yr": int(max(gaps)),
        "gap_std_yr": float(np.std(gaps)),
        "fail_year_min": fails[0],
        "fail_year_max": fails[-1],
    }


def _autocorrelation(series: np.ndarray, lag_max: int) -> list[float]:
    """Pearson autocorrelation of ``series`` at lags 1..lag_max."""
    s = series - series.mean()
    var = np.dot(s, s)
    if var == 0.0:
        return [0.0] * lag_max
    out: list[float] = []
    for lag in range(1, lag_max + 1):
        if lag >= len(s):
            out.append(0.0)
            continue
        cov = float(np.dot(s[:-lag], s[lag:]))
        out.append(cov / var)
    return out


def _best_autocorr_peak(autocorr: list[float], lag_min: int = 2) -> tuple[int, float]:
    """Best non-trivial lag + its autocorrelation value."""
    if not autocorr:
        return 0, 0.0
    best_lag = 0
    best_val = -np.inf
    for lag in range(lag_min, len(autocorr) + 1):
        v = autocorr[lag - 1]
        if v > best_val:
            best_val = v
            best_lag = lag
    return best_lag, float(best_val)


def _kernel_edge_failures(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """How many failures sit near the URA111 2099 kernel edge?"""
    edge_cutoff = KERNEL_EDGE_YEAR - KERNEL_EDGE_YEARS  # 2084
    n_fails_total = sum(1 for r in rows if not r["passes_v4_strict"])
    n_fails_edge = sum(1 for r in rows if not r["passes_v4_strict"] and r["year"] >= edge_cutoff)
    n_fails_interior = n_fails_total - n_fails_edge
    return {
        "edge_cutoff_year": edge_cutoff,
        "edge_window_years": KERNEL_EDGE_YEARS,
        "n_fails_near_edge": n_fails_edge,
        "n_fails_interior": n_fails_interior,
        "n_fails_total": n_fails_total,
    }


def _longest_pass_run(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Longest contiguous PASS stretch + its center (for the recommended
    launch_epoch in the CYCLIC / EFFECTIVELY_CYCLIC branch)."""
    best_start: int | None = None
    best_end: int | None = None
    best_len = 0
    cur_start: int | None = None
    cur_len = 0
    for r in rows:
        if r["passes_v4_strict"]:
            if cur_start is None:
                cur_start = r["year"]
            cur_len += 1
            if cur_len > best_len:
                best_len = cur_len
                best_start = cur_start
                best_end = r["year"]
        else:
            cur_start = None
            cur_len = 0
    return {
        "longest_pass_run_start_year": best_start,
        "longest_pass_run_end_year": best_end,
        "longest_pass_run_length_yr": best_len,
        "longest_pass_run_center_year": (
            (best_start + best_end) // 2 if best_start and best_end else None
        ),
    }


def _verdict(
    rows: list[dict[str, Any]],
    decade_breakdown: dict[str, dict[str, int]],
    fail_clustering: dict[str, Any],
    autocorr_pass: list[float],
    autocorr_drift: list[float],
    kernel_edge: dict[str, Any],
    longest_pass: dict[str, Any],
) -> dict[str, Any]:
    """Apply the boundary-verdict rules."""
    n_total = len(rows)
    n_pass = sum(1 for r in rows if r["passes_v4_strict"])
    pass_fraction = n_pass / n_total if n_total else 0.0

    lag_max = min(len(autocorr_pass), int(LAG_FRACTION_MAX * n_total))
    pass_lag, pass_peak = _best_autocorr_peak(autocorr_pass[:lag_max])
    drift_lag, drift_peak = _best_autocorr_peak(autocorr_drift[:lag_max])

    cyclic_pass = pass_peak >= AUTOCORR_PEAK_THRESHOLD
    cyclic_drift = drift_peak >= AUTOCORR_PEAK_THRESHOLD

    n_interior_fails = kernel_edge["n_fails_interior"]
    n_edge_fails = kernel_edge["n_fails_near_edge"]
    n_total_fails = kernel_edge["n_fails_total"]
    edge_dominated = n_total_fails > 0 and n_edge_fails / n_total_fails >= 0.75

    if cyclic_pass or cyclic_drift:
        verdict_label = "CYCLIC"
        period_yr = pass_lag if cyclic_pass else drift_lag
        physical_interpretation = (
            f"Autocorrelation peak at lag {period_yr} yr "
            f"(pass={pass_peak:.3f} drift={drift_peak:.3f}). "
            f"This is on the same order as fractional Uranus-orbital "
            f"({URANUS_ORBITAL_YEARS:.2f} yr) periods. Sub-year periodicities "
            f"(Umbriel {UMBRIEL_PERIOD_DAYS}d, Oberon {OBERON_PERIOD_DAYS}d, "
            f"U-O synodic {UMBRIEL_OBERON_SYNODIC_DAYS}d) are aliased at "
            f"annual resolution and cannot be resolved here."
        )
        next_step = (
            "CYCLIC PASS pattern detected. Fire #337 successor to compose "
            "the catalogue row as quasi_cycler with epoch_locked=true, "
            f"launch_epoch in the {longest_pass['longest_pass_run_start_year']}.."
            f"{longest_pass['longest_pass_run_end_year']} band (center "
            f"{longest_pass['longest_pass_run_center_year']}), and validity_window "
            "set to the longest contiguous PASS stretch."
        )
    elif edge_dominated and pass_fraction >= MIN_PASS_FRACTION:
        verdict_label = "EFFECTIVELY_CYCLIC"
        period_yr = None
        n_interior_years = n_total - KERNEL_EDGE_YEARS
        n_interior_pass = n_interior_years - n_interior_fails
        interior_pass_rate = n_interior_pass / n_interior_years if n_interior_years > 0 else 0.0
        physical_interpretation = (
            f"No clear sub-100-yr periodicity (autocorrelation peaks "
            f"pass={pass_peak:.3f} at lag {pass_lag} yr, drift={drift_peak:.3f} "
            f"at lag {drift_lag} yr -- both below the {AUTOCORR_PEAK_THRESHOLD} "
            f"threshold). HOWEVER, {n_edge_fails}/{n_total_fails} failures "
            f"sit in the last {KERNEL_EDGE_YEARS} yr of the URA111 kernel's "
            f"validity window. This is a kernel-edge / extrapolation artifact, "
            f"not an intrinsic resonance failure. Interior PASS rate is "
            f"{n_interior_pass}/{n_interior_years} = {interior_pass_rate:.2%}."
        )
        next_step = (
            "EFFECTIVELY CYCLIC -- failures concentrated at URA111 kernel "
            "edge, interior is uniformly PASSing. Fire #337 successor with "
            f"launch_epoch in the {longest_pass['longest_pass_run_start_year']}.."
            f"{longest_pass['longest_pass_run_end_year']} band (center "
            f"{longest_pass['longest_pass_run_center_year']}); validity_window "
            "= the longest PASS stretch; flag that beyond ~2085 a fresher "
            "Uranian satellite kernel (post-URA111) is needed."
        )
    else:
        verdict_label = "IRREGULAR"
        period_yr = None
        physical_interpretation = (
            f"No clear periodicity in either the binary PASS series "
            f"(peak {pass_peak:.3f} at lag {pass_lag} yr) or the continuous "
            f"drift_vs_v3 series (peak {drift_peak:.3f} at lag {drift_lag} yr), "
            f"both below the {AUTOCORR_PEAK_THRESHOLD} threshold. The "
            f"PASS/FAIL pattern at annual resolution looks irregular, with "
            f"{n_interior_fails} interior fails and {n_edge_fails} kernel-edge "
            f"fails. The bounded-drift signature under V4-strict is too "
            f"fragile to be mission-useful."
        )
        next_step = (
            "IRREGULAR PASS pattern. Retire to data/empty_regions.jsonl "
            "(negative-results registry) per project_negative_results_registry. "
            "Structural finding: bounded-drift signature exists under V3 "
            "(idealised) but real-ephemeris launch-epoch sensitivity breaks "
            "the resonance unpredictably; SILVER not catalogue-admissible "
            "as quasi_cycler without finer epoch resolution + a different "
            "discovery method."
        )

    return {
        "verdict_label": verdict_label,
        "pass_fraction": pass_fraction,
        "autocorr_pass_best_lag_yr": pass_lag,
        "autocorr_pass_best_peak": pass_peak,
        "autocorr_drift_best_lag_yr": drift_lag,
        "autocorr_drift_best_peak": drift_peak,
        "autocorr_threshold": AUTOCORR_PEAK_THRESHOLD,
        "lag_max_yr": lag_max,
        "kernel_edge_dominated": edge_dominated,
        "physical_interpretation": physical_interpretation,
        "next_step": next_step,
        "recommended_launch_window": {
            "start_year": longest_pass["longest_pass_run_start_year"],
            "end_year": longest_pass["longest_pass_run_end_year"],
            "center_year": longest_pass["longest_pass_run_center_year"],
            "length_yr": longest_pass["longest_pass_run_length_yr"],
        }
        if verdict_label in ("CYCLIC", "EFFECTIVELY_CYCLIC")
        else None,
    }


def main() -> int:
    sha = _git_sha()
    t0 = time.time()
    print(f"[#338-B] boundary analysis -- sha={sha}", flush=True)

    rows, _header = _read_sweep()
    print(f"[#338-B] read {len(rows)} per-year rows from {IN_JSONL.name}", flush=True)
    if not rows:
        print(f"[#338-B] FATAL: no rows in {IN_JSONL}", file=sys.stderr)
        return 1

    n_total = len(rows)
    n_pass = sum(1 for r in rows if r["passes_v4_strict"])
    n_fail = n_total - n_pass
    print(
        f"[#338-B] PASS/FAIL totals: {n_pass}/{n_fail} ({n_pass / n_total:.2%} PASS)",
        flush=True,
    )

    decades = _decade_breakdown(rows)
    print("[#338-B] per-decade breakdown:", flush=True)
    for dec in sorted(decades):
        d = decades[dec]
        print(
            f"            {dec}: PASS={d['pass']:2d} FAIL={d['fail']:2d}",
            flush=True,
        )

    fail_cluster = _fail_clustering(rows)
    print(f"[#338-B] FAIL clustering: {fail_cluster}", flush=True)

    # Autocorrelation series.
    pass_array = np.array([1.0 if r["passes_v4_strict"] else 0.0 for r in rows])
    # Use FINITE drift values for autocorrelation (replace inf with 50,000 km
    # floor + 1 so they sit just above the threshold without overwhelming
    # the variance).
    drift_array = np.array(
        [
            min(r["drift_agreement_kms_vs_v3"], 5e4 + 1.0)
            if np.isfinite(r["drift_agreement_kms_vs_v3"])
            else 5e4 + 1.0
            for r in rows
        ]
    )
    lag_max_compute = min(50, n_total - 1)
    autocorr_pass = _autocorrelation(pass_array, lag_max_compute)
    autocorr_drift = _autocorrelation(drift_array, lag_max_compute)
    pass_lag, pass_peak = _best_autocorr_peak(autocorr_pass[: int(LAG_FRACTION_MAX * n_total)])
    drift_lag, drift_peak = _best_autocorr_peak(autocorr_drift[: int(LAG_FRACTION_MAX * n_total)])
    print(
        f"[#338-B] autocorrelation pass series: best lag={pass_lag} yr, peak={pass_peak:.3f}",
        flush=True,
    )
    print(
        f"[#338-B] autocorrelation drift series: best lag={drift_lag} yr, peak={drift_peak:.3f}",
        flush=True,
    )

    kernel_edge = _kernel_edge_failures(rows)
    print(f"[#338-B] kernel-edge failures: {kernel_edge}", flush=True)

    longest_pass = _longest_pass_run(rows)
    print(f"[#338-B] longest PASS run: {longest_pass}", flush=True)

    verdict = _verdict(
        rows,
        decades,
        fail_cluster,
        autocorr_pass,
        autocorr_drift,
        kernel_edge,
        longest_pass,
    )
    print(
        f"[#338-B] VERDICT: {verdict['verdict_label']} "
        f"(pass_fraction={verdict['pass_fraction']:.2%})",
        flush=True,
    )
    print(f"[#338-B]   {verdict['physical_interpretation']}", flush=True)
    print(f"[#338-B]   next_step: {verdict['next_step']}", flush=True)

    out_rows: list[dict[str, Any]] = []
    out_rows.append(
        {
            "_meta": True,
            "task": "#338 Part B -- boundary analysis on annual V4-strict sweep",
            "successor_to_part_A": str(IN_JSONL.relative_to(ROOT)),
            "git_sha": sha,
            "n_epochs": n_total,
            "n_pass": n_pass,
            "n_fail": n_fail,
            "pass_fraction": n_pass / n_total,
            "autocorr_peak_threshold": AUTOCORR_PEAK_THRESHOLD,
            "lag_fraction_max": LAG_FRACTION_MAX,
            "kernel_edge_year": KERNEL_EDGE_YEAR,
            "kernel_edge_years_window": KERNEL_EDGE_YEARS,
            "physical_period_candidates_days_or_years": {
                "umbriel_orbital_days": UMBRIEL_PERIOD_DAYS,
                "oberon_orbital_days": OBERON_PERIOD_DAYS,
                "umbriel_oberon_synodic_days": UMBRIEL_OBERON_SYNODIC_DAYS,
                "uranus_orbital_years": URANUS_ORBITAL_YEARS,
            },
            "aliasing_note": (
                "Annual sampling cannot resolve sub-year periodicities. "
                "If the true PASS/FAIL boundary is set by Umbriel-Oberon "
                "synodic phase (5.987 days), the annual grid is aliased "
                "and a daily/weekly sweep near a known transition is "
                "required for confirmation."
            ),
        }
    )
    out_rows.append(
        {
            "kind": "decade_breakdown",
            "decades": decades,
        }
    )
    out_rows.append(
        {
            "kind": "fail_clustering",
            **fail_cluster,
        }
    )
    out_rows.append(
        {
            "kind": "autocorrelation",
            "lag_max_computed": lag_max_compute,
            "lag_max_trusted": int(LAG_FRACTION_MAX * n_total),
            "autocorr_pass_series": autocorr_pass,
            "autocorr_drift_series": autocorr_drift,
            "autocorr_pass_best_lag": pass_lag,
            "autocorr_pass_best_peak": pass_peak,
            "autocorr_drift_best_lag": drift_lag,
            "autocorr_drift_best_peak": drift_peak,
        }
    )
    out_rows.append(
        {
            "kind": "kernel_edge_breakdown",
            **kernel_edge,
        }
    )
    out_rows.append(
        {
            "kind": "longest_pass_run",
            **longest_pass,
        }
    )
    out_rows.append(
        {
            "kind": "boundary_verdict",
            **verdict,
            "elapsed_s": time.time() - t0,
        }
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for r in out_rows:
            fh.write(json.dumps(r) + "\n")

    print(f"\n[#338-B] wrote {OUT_JSONL}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
