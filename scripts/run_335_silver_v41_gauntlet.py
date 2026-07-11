"""#335 Phase 4.1 -- V4-strict gauntlet on the #327 SILVER (Umbriel-Oberon-Umbriel).

The driver
----------
:func:`cyclerfinder.data.validation.v4_uranus_strict.run_v4_uranus_strict`
(#335 Part B) re-propagates each cycle of the SILVER's Lambert tour under
scipy DOP853 + Uranus J2 + classical-moon SPICE-driven third-body
perturbations -- the full Phase 4.1 strict-real-ephemeris path.

A V4-strict PASS clears the strictest computational gate available;
catalogue admission as ``quasi_cycler`` is then unblocked (#337 = the
actual admission task; the verdict here only records the math).

Honest scope of the launch_epoch
--------------------------------
V4-strict is the FIRST V-tier verdict that is genuinely epoch-dependent
(V2/V3/V4-scipy all use circular-coplanar Kepler moons, so a single
"phase0_deg" suffices). The Uranian moon phase configuration changes the
SILVER's Lambert geometry AND the moon precession state. We therefore
run V4-strict at MULTIPLE launch epochs and report:

  * the per-epoch verdict (PASS / FAIL)
  * the worst-case agreement-vs-V3 (the strictest pass criterion)
  * a verdict ROBUSTNESS classification: ALL_PASS, ANY_FAIL, SPLIT

A SILVER that PASSes V4-strict at multiple epochs is dramatically
stronger evidence than a single-epoch PASS.

Discipline anchors
------------------
* READ-ONLY on ``data/silver_327_verified.jsonl``, ``data/silver_327_v3_verdicts.jsonl``,
  ``data/silver_327_v4_verdicts.jsonl``.
* READ-ONLY on ``src/cyclerfinder/data/validation/v[1-3]_*.py`` (Phases 1-3)
  and ``src/cyclerfinder/data/validation/v4_uranus.py`` (V4-scipy fallback).
* NO catalogue writeback. The verdict is recorded; what it means is
  documented in ``docs/notes/2026-06-16-335-v4-strict-phase41.md``.
* The verdict is whatever the math says -- PASS / FAIL -- no test-tuning.
* The verdict is recorded to ``data/silver_327_v4_strict_verdicts.jsonl``;
  the catalogue admission (#337) is a downstream task.

Run as::

    uv run python scripts/run_335_silver_v41_gauntlet.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.data.validation.v2_moontour import run_v2_moontour  # noqa: E402
from cyclerfinder.data.validation.v3_3d import V3Verdict3D, run_v3_3d  # noqa: E402
from cyclerfinder.data.validation.v4_uranus import (  # noqa: E402
    URANIAN_PERTURBER_MOONS,
    URANUS_J2,
    URANUS_R_EQ_KM,
    V4_AGREEMENT_FLOOR_KMS,
    V4_N_CYCLES_MIN,
    V4UranusVerdict,
    run_v4_uranus,
)
from cyclerfinder.data.validation.v4_uranus_strict import (  # noqa: E402
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
    V4UranusStrictVerdict,
    _verdict_to_jsonable,
    run_v4_uranus_strict,
)

SILVER_VERIFIED_JSONL = ROOT / "data" / "silver_327_verified.jsonl"
V3_VERDICTS_JSONL = ROOT / "data" / "silver_327_v3_verdicts.jsonl"
V4_SCIPY_VERDICTS_JSONL = ROOT / "data" / "silver_327_v4_verdicts.jsonl"
OUT_JSONL = ROOT / "data" / "silver_327_v4_strict_verdicts.jsonl"

# Stored SILVER row fields -- same constants as scripts/run_332_*.py
# (sourced from data/silver_327_verified.jsonl + silver_327_v3_verdicts.jsonl,
# READ-ONLY in this script -- inputs only).
SILVER_ID = "repeated-moon-uranus-00000041"
SILVER_SEQ: tuple[str, ...] = ("Umbriel", "Oberon", "Umbriel")
SILVER_VINF: tuple[float, ...] = (
    0.9199258810725036,
    0.9604309791298091,
    0.8946936085078939,
)
SILVER_TOF: tuple[float, ...] = (14.940560615336594, 14.940560615336594)
SILVER_REL_OFF_DEG = 180.0
SILVER_NREV: tuple[int, ...] = (1, 1)
SILVER_PHASE0_DEG = 29.999999999999996

# Launch epoch sweep: the SILVER's V2/V3/V4-scipy work was epoch-blind, so
# there's no SOURCED "right" epoch to pick. We sample 3 epochs spanning the
# ura111 coverage to test launch-epoch sensitivity:
#   - 2000-01-15: J2000 reference (matches the smoke test in tests/)
#   - 2030-06-21: 30 years later (different Uranian satellite configuration)
#   - 2050-12-31: 50 years later (near the kernel's 2099 endpoint)
# A PASS robust across these 3 is dramatically stronger than single-epoch.
LAUNCH_EPOCHS_UTC: tuple[str, ...] = (
    "2000-01-15T00:00:00",
    "2030-06-21T00:00:00",
    "2050-12-31T00:00:00",
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _build_v3_v4scipy(
    n_cycles: int,
) -> tuple[V3Verdict3D, V4UranusVerdict]:
    """V2 -> V3 -> V4-scipy chain (epoch-blind, same as #331/#332)."""
    v2 = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        n_cycles=n_cycles,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
        notes=f"#335 V4-strict input chain, n_cycles={n_cycles}",
    )
    v3 = run_v3_3d(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        v2_verdict=v2,
        n_cycles=n_cycles,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
        notes=f"#335 V4-strict input chain, n_cycles={n_cycles}",
    )
    v4 = run_v4_uranus(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        v3_verdict=v3,
        n_cycles=n_cycles,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
        notes=f"#335 V4-strict input chain, n_cycles={n_cycles}",
    )
    return v3, v4


def _epoch_safe_label(epoch: str) -> str:
    """Filesystem-safe epoch label (strip colons + T)."""
    return epoch.replace(":", "").replace("T", "-")


def main() -> int:
    sha = _git_sha()
    t0 = time.time()
    print(f"[#335] V4-strict gauntlet Phase 4.1 -- sha={sha}", flush=True)
    print(f"[#335] candidate = {SILVER_ID}", flush=True)
    print(f"[#335] sequence  = {SILVER_SEQ}, n_rev = {SILVER_NREV}", flush=True)
    print(
        f"[#335] SPICE kernels = LSK + PCK + ura111.bsp ({DEFAULT_URA_PATH})",
        flush=True,
    )
    print(
        f"[#335] launch epochs = {LAUNCH_EPOCHS_UTC} (testing epoch sensitivity)",
        flush=True,
    )

    # Sanity-check kernels are present BEFORE running the long chain.
    for p in (DEFAULT_LSK_PATH, DEFAULT_PCK_PATH, DEFAULT_URA_PATH):
        if not p.exists():
            print(f"[#335] FATAL: SPICE kernel missing: {p}", file=sys.stderr)
            return 1

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#335 V4-strict Phase 4.1 gauntlet -- Umbriel-Oberon-Umbriel SILVER",
            "candidate_id": SILVER_ID,
            "source_jsonl": str(SILVER_VERIFIED_JSONL.relative_to(ROOT)),
            "v3_verdicts_jsonl": str(V3_VERDICTS_JSONL.relative_to(ROOT)),
            "v4_scipy_verdicts_jsonl": str(V4_SCIPY_VERDICTS_JSONL.relative_to(ROOT)),
            "stored_silver": {
                "sequence": list(SILVER_SEQ),
                "vinf_per_encounter_kms": list(SILVER_VINF),
                "tof_days": list(SILVER_TOF),
                "rel_offset_deg": SILVER_REL_OFF_DEG,
                "phase0_deg": SILVER_PHASE0_DEG,
                "n_rev": list(SILVER_NREV),
            },
            "v4_strict_model": {
                "integrator": "scipy DOP853",
                "uranus_j2": URANUS_J2,
                "uranus_r_eq_km": URANUS_R_EQ_KM,
                "perturber_moons": list(URANIAN_PERTURBER_MOONS),
                "j2_source": "Jacobson 2014, AJ 148:76, Table 4",
                "spice_kernels_used": [
                    str(DEFAULT_LSK_PATH),
                    str(DEFAULT_PCK_PATH),
                    str(DEFAULT_URA_PATH),
                ],
                "spice_kernel_source": (
                    "JPL/NAIF generic_kernels/spk/satellites/a_old_versions/ura111.bsp; "
                    "5 classical Uranian moons (Miranda/Ariel/Umbriel/Titania/Oberon) "
                    "+ Uranus + Earth + Sun, 1900-2099 ET"
                ),
            },
            "driver_floors": {
                "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                "n_cycles_min": V4_N_CYCLES_MIN,
            },
            "launch_epochs_utc": list(LAUNCH_EPOCHS_UTC),
            "n_cycles_grid": [3, 5, 10],
            "git_sha": sha,
        }
    )

    n_cycles_grid = (3, 5, 10)
    # Cache the V3 + V4-scipy chains per n_cycles (epoch-blind so once per nc).
    chains: dict[int, tuple[V3Verdict3D, V4UranusVerdict]] = {}
    for nc in n_cycles_grid:
        print(f"[#335] V2->V3->V4-scipy chain at n_cycles={nc}...", flush=True)
        t_chain = time.time()
        chains[nc] = _build_v3_v4scipy(nc)
        print(f"  -> chain ready (elapsed {time.time() - t_chain:.1f}s)", flush=True)

    # Headline per-epoch verdicts (n_cycles=3, 5, 10).
    headline_v4_strict: dict[tuple[str, int], V4UranusStrictVerdict] = {}

    for epoch in LAUNCH_EPOCHS_UTC:
        print(f"\n[#335] === LAUNCH EPOCH {epoch} ===", flush=True)
        for nc in n_cycles_grid:
            v3, v4_scipy = chains[nc]
            print(f"[#335] run_v4_uranus_strict(epoch={epoch}, n_cycles={nc})...", flush=True)
            t_run = time.time()
            v4s = run_v4_uranus_strict(
                SILVER_ID,
                SILVER_SEQ,
                SILVER_VINF,
                SILVER_TOF,
                SILVER_REL_OFF_DEG,
                epoch,
                None,
                v3_verdict=v3,
                v4_scipy_verdict=v4_scipy,
                n_cycles=nc,
                n_revs=SILVER_NREV,
                notes=f"#335 phase-4.1 V4-strict scan, epoch={epoch}, n_cycles={nc}",
            )
            headline_v4_strict[(epoch, nc)] = v4s
            print(
                f"  -> passes_v4_strict={v4s.passes_v4_strict} | "
                f"bounded={v4s.bounded_drift_survives} | "
                f"completed={v4s.n_cycles_propagated}/{nc} | "
                f"agreement_vs_v3={v4s.drift_agreement_kms_vs_v3:.3e} km | "
                f"agreement_vs_v4_scipy={v4s.drift_agreement_kms_vs_v4_scipy:.3e} km | "
                f"e_body1={v4s.eccentricity_used_e_body1:.5f} | "
                f"elapsed={time.time() - t_run:.1f}s",
                flush=True,
            )
            for c in v4s.per_cycle:
                print(
                    f"     cycle {c.cycle_index}: "
                    f"V4_strict={c.rendezvous_drift_kms_v4_strict:.3e} | "
                    f"V4_scipy={c.rendezvous_drift_kms_v4_scipy:.3e} | "
                    f"V3={c.rendezvous_drift_kms_v3:.3e} | "
                    f"|V4s-V3|={c.agreement_kms_vs_v3:.3e} | "
                    f"|V4s-V4|={c.agreement_kms_vs_v4_scipy:.3e}",
                    flush=True,
                )
            row = _verdict_to_jsonable(v4s)
            row["epoch_label"] = _epoch_safe_label(epoch)
            rows.append(row)

    # Per-epoch + cross-epoch summary.
    epoch_pass_grid: dict[str, dict[int, bool]] = {
        epoch: {nc: headline_v4_strict[(epoch, nc)].passes_v4_strict for nc in n_cycles_grid}
        for epoch in LAUNCH_EPOCHS_UTC
    }
    all_pass = all(
        headline_v4_strict[(epoch, nc)].passes_v4_strict
        for epoch in LAUNCH_EPOCHS_UTC
        for nc in n_cycles_grid
    )
    any_pass = any(
        headline_v4_strict[(epoch, nc)].passes_v4_strict
        for epoch in LAUNCH_EPOCHS_UTC
        for nc in n_cycles_grid
    )
    if all_pass:
        verdict_label = "PASS"
        robustness = "ROBUST_ACROSS_EPOCHS"
        next_step = (
            "10 of 10 computational gates cleared (closure, DOP853 cross-check, "
            "physical-sanity, lit-novelty, ML flagger, V1 3D, V2 moontour, "
            "V3 IAS15, V4-scipy J2+nbody fallback, V4-strict full-SPICE). "
            "The SILVER has cleared the entire computational gauntlet at "
            "multiple launch epochs -- the bounded-drift signature is robust "
            "to real Uranian satellite eccentricity, inclination, AND secular "
            "precession across decade-scale epoch variation. Catalogue "
            "admission as quasi_cycler is unblocked: fire #337 to compose "
            "the catalogue row with V0 sourced from the SILVER's full "
            "provenance chain (V1-V4-strict verdicts cited; validity_window "
            "spans the launch_epoch grid)."
        )
    elif any_pass:
        verdict_label = "MIXED"
        robustness = "EPOCH_DEPENDENT"
        next_step = (
            "Some launch epochs PASS V4-strict and others FAIL. The SILVER's "
            "bounded-drift signature is epoch-dependent on real Uranian "
            "satellite phasing -- expected for a SILVER (epoch-blind discovery) "
            "that hasn't been epoch-locked to a specific Uranian-moon "
            "configuration. Determine the PASS sub-window (validity_window) "
            "before catalogue admission. The signature is real; the geometry "
            "of when it works needs more characterisation."
        )
    else:
        verdict_label = "FAIL"
        robustness = "FAILS_ALL_EPOCHS"
        next_step = (
            "Bounded-drift signature does NOT survive real Uranian satellite "
            "eccentricity / inclination / secular precession at ANY of the "
            "tested launch epochs. Retire to negative-results registry (#172) "
            "recording the V4-strict perturbation order that breaks the "
            "signature. The V3 + V4-scipy bounded-drift was a model-class "
            "artifact that the SPICE-ephemeris reality kills."
        )

    # Per-epoch worst-case agreement (for the validity_window discussion).
    worst_agreement_vs_v3: dict[str, float] = {
        epoch: max(
            headline_v4_strict[(epoch, nc)].drift_agreement_kms_vs_v3 for nc in n_cycles_grid
        )
        for epoch in LAUNCH_EPOCHS_UTC
    }
    worst_agreement_vs_v4_scipy: dict[str, float] = {
        epoch: max(
            headline_v4_strict[(epoch, nc)].drift_agreement_kms_vs_v4_scipy for nc in n_cycles_grid
        )
        for epoch in LAUNCH_EPOCHS_UTC
    }

    rows.append(
        {
            "_meta": True,
            "kind": "headline",
            "candidate_id": SILVER_ID,
            "verdict_label": verdict_label,
            "robustness": robustness,
            "epoch_pass_grid": {
                epoch: {str(nc): epoch_pass_grid[epoch][nc] for nc in n_cycles_grid}
                for epoch in LAUNCH_EPOCHS_UTC
            },
            "worst_agreement_kms_vs_v3_per_epoch": worst_agreement_vs_v3,
            "worst_agreement_kms_vs_v4_scipy_per_epoch": worst_agreement_vs_v4_scipy,
            "writeback_to_catalogue": False,
            "next_step": next_step,
            "elapsed_s": time.time() - t0,
        }
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")

    print(f"\n[#335] wrote {OUT_JSONL}", flush=True)
    print(
        f"[#335] HEADLINE verdict: {verdict_label} (robustness={robustness}, "
        f"elapsed {time.time() - t0:.1f}s)",
        flush=True,
    )
    print(f"[#335] next step: {next_step}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
