"""#566 -- V2->V3->V4->V4-strict gauntlet on 5 representative #563 closures.

Background
----------
Task #312 (catalogue row ``umbriel-oberon-1-1-uranian-quasi-cycler-2026``) is
a confirmed, literature-novel Uranian quasi_cycler. #558-#565 established
that #312 is 1 of 30 exact symmetric periodic closures spread across the 6
non-Miranda Uranian moon-pair directions enumerated in
``data/enumerate_563_symmetric_closures.jsonl``, and recommended running the
same computational gauntlet #312 itself cleared (#330 V2 -> #331 V3 -> #332
V4 -> #335 V4-strict) on 5 more representatives, one per other pair
direction.

The driver
----------
The underlying gauntlet functions are fully generic over
``sequence``/``vinf_tuple_kms``/``leg_tofs_days``/``rel_offset_deg``/
``n_revs``/``phase0_deg``:

* :func:`cyclerfinder.data.validation.v2_moontour.run_v2_moontour`
* :func:`cyclerfinder.data.validation.v3_3d.run_v3_3d`
* :func:`cyclerfinder.data.validation.v4_uranus.run_v4_uranus`
* :func:`cyclerfinder.data.validation.v4_uranus_strict.run_v4_uranus_strict`

This script does NOT reimplement any of them -- it imports and calls the same
functions #330-#338 already used, looped over 5 new candidate parameter sets
instead of the #327 Umbriel-Oberon SILVER. ``phase0_deg`` reuses the fixed
constant already used across #330/#335/#338 (rotation-redundant per #558).

Candidate provenance (verified against the jsonl directly, NOT transcribed
blind from the #566 OUTSTANDING bullet -- the bullet's "Row N" labels do not
match ``enumerate_563_symmetric_closures.jsonl`` line numbers, but every
numeric field -- sequence, tof_days, vinf, rel_offset_deg, n_rev,
n_commensurate_int -- was independently matched to the real jsonl row it came
from; see ``source_jsonl_line`` on each candidate below):

1. jsonl line 57 (Titania-Oberon-Titania) -- MANDATORY, carries the
   literature-clearance obligation from #565.
2. jsonl line 2  (Ariel-Umbriel-Ariel)
3. jsonl line 12 (Ariel-Titania-Ariel)
4. jsonl line 18 (Ariel-Oberon-Ariel)
5. jsonl line 26 (Umbriel-Titania-Umbriel)

V4-strict epoch choice
-----------------------
V4-strict is the first epoch-dependent tier. Per #566's scope, this script
runs ONE representative launch epoch per candidate (NOT #338's full annual
sweep, which is explicitly out of scope). We use 2000-06-21T00:00:00 -- the
same reference epoch #338 originally anchored on for #312 and that #559's
731-epoch daily scan later corroborated as a PASS point in both 2000 and
2030 (see OUTSTANDING #559). Reusing an already-established reference epoch
avoids introducing a fresh, unmotivated epoch choice for this pass.

Known caveat carried into the results (NOT a capability blocker)
------------------------------------------------------------------
``run_v4_uranus_strict``'s audit-only fields ``eccentricity_used_e_umbriel``
/ ``eccentricity_used_e_oberon`` / ``inclination_used_deg_umbriel`` /
``inclination_used_deg_oberon`` are hardcoded in
``src/cyclerfinder/data/validation/v4_uranus_strict.py`` (lines ~575-577) to
always SPICE-sample Umbriel and Oberon specifically, regardless of the
candidate's actual ``sequence``. The core physics propagation
(``_cycle_v4_strict``) is fully generic over ``sequence`` and
``perturber_moons`` (uses ``_moon_state_spice(moon, ...)`` per actual moon
name) -- so the PASS/FAIL verdict itself is unaffected. Only these 4
descriptive audit fields are mislabeled/irrelevant for the 4 non-Umbriel-
Oberon candidates in this run; this script records them as-is (for
traceability) but the per-stage numeric detail below should be read via the
generic drift/agreement fields, not the moon-specific eccentricity fields,
for any candidate other than a genuine Umbriel-Oberon pair.

Discipline anchors
-------------------
* READ-ONLY on ``data/enumerate_563_symmetric_closures.jsonl``.
* READ-ONLY on ``src/cyclerfinder/data/validation/v[2-4]*.py`` -- no new
  capability, no driver modification.
* NO catalogue writeback -- ``data/catalogue.yaml`` is untouched. That
  remains a deliberate, separate step regardless of how these 5 candidates
  fare.
* NOT the #338-style annual DOY sweep -- single representative epoch only.

Run as::

    uv run python scripts/run_566_gauntlet_five_representatives.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.data.validation.v2_moontour import (  # noqa: E402
    V2_MOONTOUR_CLOSURE_FLOOR_KMS,
    V2_MOONTOUR_DRIFT_FLOOR_KMS,
    V2_MOONTOUR_N_CYCLES_MIN,
    V2MoontourVerdict,
    run_v2_moontour,
)
from cyclerfinder.data.validation.v3_3d import (  # noqa: E402
    V3_AGREEMENT_FLOOR_KMS,
    V3_N_CYCLES_MIN,
    V3Verdict3D,
    run_v3_3d,
)
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

SOURCE_JSONL = ROOT / "data" / "enumerate_563_symmetric_closures.jsonl"
OUT_JSONL = ROOT / "data" / "gauntlet_566_five_representatives.jsonl"

# #558-established rotation-redundant constant, reused verbatim from
# #330/#335/#338 -- not re-derived here.
PHASE0_DEG = 29.999999999999996

# Single representative launch epoch for V4-strict (see module docstring for
# rationale -- reuses the #338/#559-established reference epoch).
LAUNCH_EPOCH_UTC = "2000-06-21T00:00:00"

N_CYCLES_GRID: tuple[int, ...] = (3, 5, 10)


@dataclass(frozen=True)
class Candidate:
    """One #563 candidate closure, verified against the source jsonl."""

    candidate_id: str
    label: str
    mandatory: bool
    source_jsonl_line: int
    sequence: tuple[str, ...]
    vinf_kms: tuple[float, ...]
    tof_days: tuple[float, ...]
    rel_offset_deg: float
    n_revs: tuple[int, ...]
    n_commensurate_int: int
    t_syn_days: float
    source_residual_kms: float


CANDIDATES: tuple[Candidate, ...] = (
    Candidate(
        candidate_id="enum563-line57-titania-oberon-titania",
        label="Titania-Oberon-Titania (MANDATORY -- carries literature-clearance obligation)",
        mandatory=True,
        source_jsonl_line=57,
        sequence=("Titania", "Oberon", "Titania"),
        vinf_kms=(2.161767816675378, 1.9680495724521725, 2.1617678166753755),
        tof_days=(12.316046445872583, 12.316046445872583),
        rel_offset_deg=180.0,
        n_revs=(0, 0),
        n_commensurate_int=1,
        t_syn_days=24.632092891745167,
        source_residual_kms=2.6645352591003757e-15,
    ),
    Candidate(
        candidate_id="enum563-line2-ariel-umbriel-ariel",
        label="Ariel-Umbriel-Ariel",
        mandatory=False,
        source_jsonl_line=2,
        sequence=("Ariel", "Umbriel", "Ariel"),
        vinf_kms=(0.979040480994661, 1.3004234339628207, 0.9790404809946573),
        tof_days=(3.216088179066208, 3.216088179066208),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
        n_commensurate_int=1,
        t_syn_days=6.432176358132416,
        source_residual_kms=1.1324274851176597e-14,
    ),
    Candidate(
        candidate_id="enum563-line12-ariel-titania-ariel",
        label="Ariel-Titania-Ariel",
        mandatory=False,
        source_jsonl_line=12,
        sequence=("Ariel", "Titania", "Ariel"),
        vinf_kms=(1.2306411593828481, 1.7185773183747601, 1.2306411593828457),
        tof_days=(5.320895317317783, 5.320895317317783),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
        n_commensurate_int=3,
        t_syn_days=3.547263544878522,
        source_residual_kms=2.4424906541753444e-15,
    ),
    Candidate(
        candidate_id="enum563-line18-ariel-oberon-ariel",
        label="Ariel-Oberon-Ariel",
        mandatory=False,
        source_jsonl_line=18,
        sequence=("Ariel", "Oberon", "Ariel"),
        vinf_kms=(1.520866047614147, 1.8285940380726622, 1.5208660476141462),
        tof_days=(7.751820498940574, 7.751820498940574),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
        n_commensurate_int=5,
        t_syn_days=3.1007281995762295,
        source_residual_kms=1.9984014443252818e-15,
    ),
    Candidate(
        candidate_id="enum563-line26-umbriel-titania-umbriel",
        label="Umbriel-Titania-Umbriel",
        mandatory=False,
        source_jsonl_line=26,
        sequence=("Umbriel", "Titania", "Umbriel"),
        vinf_kms=(1.2295656768416439, 1.0058255988095806, 1.229565676841644),
        tof_days=(3.9544738760575804, 3.9544738760575804),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
        n_commensurate_int=1,
        t_syn_days=7.908947752115161,
        source_residual_kms=4.440892098500626e-15,
    ),
)

_REGION_ID = "uranian-563-five-representative-symmetric-closures-2026-07-11"
_METHOD = MethodCapability(
    genome=(
        "5 representative exact symmetric periodic closures (one per non-Miranda "
        "Uranian moon-pair direction other than Umbriel-Oberon) from the #563 "
        "enumeration, run through the existing frozen V2->V3->V4->V4-strict "
        "validation chain -- read-only validation reporting, no genome/corrector change"
    ),
    corrector=(
        "existing V2->V3->V4-scipy->V4-strict chain "
        "(validation/v2_moontour.py, v3_3d.py, v4_uranus.py, v4_uranus_strict.py), unmodified"
    ),
    capability_tags=frozenset(
        {"cr3bp", "real-ephemeris", "v4-strict", "uranian", "validation-reporting", "moontour"}
    ),
    git_sha="working-tree",
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _v2_row(candidate_id: str, verdict: V2MoontourVerdict) -> dict[str, Any]:
    return {
        "kind": "moontour_v2_verdict",
        "candidate_id": candidate_id,
        "sequence": list(verdict.sequence),
        "n_cycles_requested": verdict.n_cycles_requested,
        "n_cycles_completed": verdict.n_cycles_completed,
        "drift_floor_kms": verdict.drift_floor_kms,
        "closure_floor_kms": verdict.closure_floor_kms,
        "n_cycles_min": verdict.n_cycles_min,
        "max_drift_kms": verdict.max_drift_kms,
        "max_drift_seconds": verdict.max_drift_seconds,
        "max_closure_residual_kms": verdict.max_closure_residual_kms,
        "passes_v2": verdict.passes_v2,
        "notes": verdict.notes,
    }


def _v3_row(candidate_id: str, verdict: V3Verdict3D) -> dict[str, Any]:
    return {
        "kind": "moontour_v3_verdict",
        "candidate_id": candidate_id,
        "sequence": list(verdict.sequence),
        "n_cycles_propagated": verdict.n_cycles_propagated,
        "integrator": verdict.integrator,
        "drift_agreement_kms": verdict.drift_agreement_kms,
        "v3_v2_agreement_floor_kms": verdict.v3_v2_agreement_floor_kms,
        "passes_v3": verdict.passes_v3,
        "notes": verdict.notes,
    }


def _v4_row(candidate_id: str, verdict: V4UranusVerdict) -> dict[str, Any]:
    return {
        "kind": "moontour_v4_verdict",
        "candidate_id": candidate_id,
        "sequence": list(verdict.sequence),
        "n_cycles_propagated": verdict.n_cycles_propagated,
        "integrator": verdict.integrator,
        "drift_agreement_kms": verdict.drift_agreement_kms,
        "v4_v3_agreement_floor_kms": verdict.v4_v3_agreement_floor_kms,
        "bounded_drift_survives": verdict.bounded_drift_survives,
        "passes_v4": verdict.passes_v4,
        "notes": verdict.notes,
    }


def _run_chain_at_ncycles(
    cand: Candidate, nc: int
) -> tuple[V2MoontourVerdict, V3Verdict3D, V4UranusVerdict]:
    v2 = run_v2_moontour(
        cand.candidate_id,
        cand.sequence,
        cand.vinf_kms,
        cand.tof_days,
        cand.rel_offset_deg,
        None,
        n_cycles=nc,
        n_revs=cand.n_revs,
        phase0_deg=PHASE0_DEG,
        notes=f"#566 gauntlet, {cand.label}, n_cycles={nc}",
    )
    v3 = run_v3_3d(
        cand.candidate_id,
        cand.sequence,
        cand.vinf_kms,
        cand.tof_days,
        cand.rel_offset_deg,
        None,
        v2_verdict=v2,
        n_cycles=nc,
        n_revs=cand.n_revs,
        phase0_deg=PHASE0_DEG,
        notes=f"#566 gauntlet, {cand.label}, n_cycles={nc}",
    )
    v4 = run_v4_uranus(
        cand.candidate_id,
        cand.sequence,
        cand.vinf_kms,
        cand.tof_days,
        cand.rel_offset_deg,
        None,
        v3_verdict=v3,
        n_cycles=nc,
        n_revs=cand.n_revs,
        phase0_deg=PHASE0_DEG,
        notes=f"#566 gauntlet, {cand.label}, n_cycles={nc}",
    )
    return v2, v3, v4


def main() -> int:
    preflight_search(
        task_no=566,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(CANDIDATES),
        override_reason=(
            "read-only validation-reporting gauntlet on 5 already-enumerated #563 "
            "candidate closures using the existing frozen V2->V3->V4->V4-strict "
            "pipeline unmodified -- not a discovery sweep; 5 candidates x 4 gauntlet "
            "stages is well under the timing-pilot threshold and the analogous #338 "
            "100-epoch single-SILVER sweep already ran in 42s, so this is seconds-scale."
        ),
    )

    sha = _git_sha()
    t0 = time.time()
    print(f"[#566] gauntlet on 5 #563 representative closures -- sha={sha}", flush=True)
    print(f"[#566] phase0_deg={PHASE0_DEG} (fixed, #558 rotation-redundant)", flush=True)
    print(
        f"[#566] V4-strict launch epoch = {LAUNCH_EPOCH_UTC} (single representative epoch)",
        flush=True,
    )

    for p in (DEFAULT_LSK_PATH, DEFAULT_PCK_PATH, DEFAULT_URA_PATH):
        if not p.exists():
            print(f"[#566] FATAL: SPICE kernel missing: {p}", file=sys.stderr)
            return 1

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#566 gauntlet on 5 #563 representative symmetric closures",
            "source_jsonl": str(SOURCE_JSONL.relative_to(ROOT)),
            "phase0_deg": PHASE0_DEG,
            "launch_epoch_utc": LAUNCH_EPOCH_UTC,
            "n_cycles_grid": list(N_CYCLES_GRID),
            "driver_floors": {
                "v2_drift_floor_kms": V2_MOONTOUR_DRIFT_FLOOR_KMS,
                "v2_closure_floor_kms": V2_MOONTOUR_CLOSURE_FLOOR_KMS,
                "v2_n_cycles_min": V2_MOONTOUR_N_CYCLES_MIN,
                "v3_agreement_floor_kms": V3_AGREEMENT_FLOOR_KMS,
                "v3_n_cycles_min": V3_N_CYCLES_MIN,
                "v4_agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                "v4_n_cycles_min": V4_N_CYCLES_MIN,
            },
            "v4_strict_model": {
                "integrator": "scipy DOP853",
                "uranus_j2": URANUS_J2,
                "uranus_r_eq_km": URANUS_R_EQ_KM,
                "perturber_moons": list(URANIAN_PERTURBER_MOONS),
                "spice_kernels_used": [
                    str(DEFAULT_LSK_PATH),
                    str(DEFAULT_PCK_PATH),
                    str(DEFAULT_URA_PATH),
                ],
            },
            "known_caveat": (
                "run_v4_uranus_strict's audit-only eccentricity/inclination fields "
                "are hardcoded to sample Umbriel/Oberon regardless of the candidate's "
                "actual sequence; physics propagation is fully generic over sequence "
                "(see module docstring). Does not affect passes_v4_strict."
            ),
            "candidates": [
                {
                    "candidate_id": c.candidate_id,
                    "label": c.label,
                    "mandatory": c.mandatory,
                    "source_jsonl_line": c.source_jsonl_line,
                    "sequence": list(c.sequence),
                    "vinf_kms": list(c.vinf_kms),
                    "tof_days": list(c.tof_days),
                    "rel_offset_deg": c.rel_offset_deg,
                    "n_revs": list(c.n_revs),
                    "n_commensurate_int": c.n_commensurate_int,
                    "t_syn_days": c.t_syn_days,
                    "source_residual_kms": c.source_residual_kms,
                }
                for c in CANDIDATES
            ],
            "git_sha": sha,
        }
    )

    summary: list[dict[str, Any]] = []

    for cand in CANDIDATES:
        t_cand = time.time()
        print(f"\n[#566] === candidate {cand.candidate_id} ({cand.label}) ===", flush=True)
        print(
            f"[#566] sequence={cand.sequence}, n_rev={cand.n_revs}, "
            f"rel_off={cand.rel_offset_deg}deg, vinf={cand.vinf_kms}, tof={cand.tof_days}",
            flush=True,
        )

        chains: dict[int, tuple[V2MoontourVerdict, V3Verdict3D, V4UranusVerdict]] = {}
        for nc in N_CYCLES_GRID:
            print(f"[#566]   V2->V3->V4 chain at n_cycles={nc}...", flush=True)
            v2, v3, v4 = _run_chain_at_ncycles(cand, nc)
            chains[nc] = (v2, v3, v4)
            print(
                f"     v2: passes={v2.passes_v2} drift={v2.max_drift_kms:.3e}km "
                f"closure={v2.max_closure_residual_kms:.3e}km/s",
                flush=True,
            )
            print(
                f"     v3: passes={v3.passes_v3} agreement={v3.drift_agreement_kms:.3e}km",
                flush=True,
            )
            print(
                f"     v4: passes={v4.passes_v4} bounded={v4.bounded_drift_survives} "
                f"agreement={v4.drift_agreement_kms:.3e}km",
                flush=True,
            )
            rows.append(_v2_row(cand.candidate_id, v2))
            rows.append(_v3_row(cand.candidate_id, v3))
            rows.append(_v4_row(cand.candidate_id, v4))

        v2_all_pass = all(chains[nc][0].passes_v2 for nc in N_CYCLES_GRID)
        v3_all_pass = all(chains[nc][1].passes_v3 for nc in N_CYCLES_GRID)
        v4_all_pass = all(chains[nc][2].passes_v4 for nc in N_CYCLES_GRID)
        # Mirror #330's quasi_cycler_3cycle test verbatim (the same test that
        # let the #312/#327 SILVER itself proceed past a strict V2 drift-
        # floor FAIL to V3/V4/V4-strict and eventual quasi_cycler catalogue
        # admission): every cycle's Lambert converged AND the per-cycle
        # V_inf-continuity closure residual stayed comfortably bounded
        # (< 0.5 km/s) even though the 50,000 km drift floor was blown.
        v2_quasi_bounded = (not v2_all_pass) and all(
            chains[nc][0].max_closure_residual_kms < 0.5 for nc in N_CYCLES_GRID
        )

        # V4-strict: single representative epoch, n_cycles grid.
        v4_strict_by_nc: dict[int, V4UranusStrictVerdict] = {}
        print(f"[#566]   V4-strict at epoch={LAUNCH_EPOCH_UTC}...", flush=True)
        for nc in N_CYCLES_GRID:
            _, v3, v4 = chains[nc]
            v4s = run_v4_uranus_strict(
                cand.candidate_id,
                cand.sequence,
                cand.vinf_kms,
                cand.tof_days,
                cand.rel_offset_deg,
                LAUNCH_EPOCH_UTC,
                None,
                v3_verdict=v3,
                v4_scipy_verdict=v4,
                n_cycles=nc,
                n_revs=cand.n_revs,
                notes=f"#566 gauntlet, {cand.label}, epoch={LAUNCH_EPOCH_UTC}, n_cycles={nc}",
            )
            v4_strict_by_nc[nc] = v4s
            print(
                f"     v4_strict(n_cycles={nc}): passes={v4s.passes_v4_strict} "
                f"bounded={v4s.bounded_drift_survives} "
                f"agreement_vs_v3={v4s.drift_agreement_kms_vs_v3:.3e}km "
                f"agreement_vs_v4_scipy={v4s.drift_agreement_kms_vs_v4_scipy:.3e}km",
                flush=True,
            )
            row = _verdict_to_jsonable(v4s)
            row["epoch_label"] = LAUNCH_EPOCH_UTC.replace(":", "").replace("T", "-")
            rows.append(row)

        v4_strict_all_pass = all(v4_strict_by_nc[nc].passes_v4_strict for nc in N_CYCLES_GRID)

        # V2 status label mirrors #330's own three-way verdict_label
        # (PASS / FAIL_QUASI_BOUNDED / FAIL) rather than treating the strict
        # 50,000 km drift floor as a hard gate -- the quasi_cycler class
        # (which #312 itself belongs to) is BY DEFINITION admitted through
        # the FAIL_QUASI_BOUNDED path, not a strict V2 PASS.
        if v2_all_pass:
            v2_status = "PASS"
        elif v2_quasi_bounded:
            v2_status = "FAIL_QUASI_BOUNDED"
        else:
            v2_status = "FAIL_UNBOUNDED"

        # V3/V4/V4-strict are independently computed regardless of the V2
        # status (matching #331/#332/#335's own unconditional-chain
        # behaviour) -- they test agreement with V2's ACTUAL trajectory, not
        # whether V2 cleared its own drift floor. "Proceeding" past V2 to
        # these stages is only substantively meaningful when V2 is PASS or
        # FAIL_QUASI_BOUNDED (a true FAIL_UNBOUNDED V2 makes the downstream
        # agreement checks moot even if they numerically pass, since V2's
        # own trajectory is no longer physically bounded).
        v2_admits_downstream = v2_status in ("PASS", "FAIL_QUASI_BOUNDED")

        if v2_admits_downstream and v3_all_pass and v4_all_pass and v4_strict_all_pass:
            highest_stage = "V4_STRICT"
            chain_verdict = "PASS_AS_QUASI_CYCLER" if v2_status == "FAIL_QUASI_BOUNDED" else "PASS"
        elif v2_admits_downstream and v3_all_pass and v4_all_pass:
            highest_stage = "V4"
            chain_verdict = "FAIL_AT_V4_STRICT"
        elif v2_admits_downstream and v3_all_pass:
            highest_stage = "V3"
            chain_verdict = "FAIL_AT_V4"
        elif v2_admits_downstream:
            highest_stage = "V2"
            chain_verdict = "FAIL_AT_V3"
        else:
            highest_stage = "NONE"
            chain_verdict = "FAIL_AT_V2_UNBOUNDED"

        cand_summary = {
            "candidate_id": cand.candidate_id,
            "label": cand.label,
            "mandatory": cand.mandatory,
            "v2_all_pass": v2_all_pass,
            "v2_quasi_bounded": v2_quasi_bounded,
            "v2_status": v2_status,
            "v2_per_nc": {
                str(nc): {
                    "passes_v2": chains[nc][0].passes_v2,
                    "max_drift_kms": chains[nc][0].max_drift_kms,
                    "max_closure_residual_kms": chains[nc][0].max_closure_residual_kms,
                }
                for nc in N_CYCLES_GRID
            },
            "v3_all_pass": v3_all_pass,
            "v3_per_nc": {
                str(nc): {
                    "passes_v3": chains[nc][1].passes_v3,
                    "drift_agreement_kms": chains[nc][1].drift_agreement_kms,
                }
                for nc in N_CYCLES_GRID
            },
            "v4_all_pass": v4_all_pass,
            "v4_per_nc": {
                str(nc): {
                    "passes_v4": chains[nc][2].passes_v4,
                    "bounded_drift_survives": chains[nc][2].bounded_drift_survives,
                    "drift_agreement_kms": chains[nc][2].drift_agreement_kms,
                }
                for nc in N_CYCLES_GRID
            },
            "v4_strict_all_pass": v4_strict_all_pass,
            "v4_strict_per_nc": {
                str(nc): {
                    "passes_v4_strict": v4_strict_by_nc[nc].passes_v4_strict,
                    "bounded_drift_survives": v4_strict_by_nc[nc].bounded_drift_survives,
                    "agreement_kms_vs_v3": v4_strict_by_nc[nc].drift_agreement_kms_vs_v3,
                    "agreement_kms_vs_v4_scipy": v4_strict_by_nc[
                        nc
                    ].drift_agreement_kms_vs_v4_scipy,
                }
                for nc in N_CYCLES_GRID
            },
            "highest_stage_all_pass": highest_stage,
            "chain_verdict": chain_verdict,
            "elapsed_s": time.time() - t_cand,
        }
        summary.append(cand_summary)
        print(
            f"[#566] candidate {cand.candidate_id} -> {chain_verdict} "
            f"(highest stage all-pass: {highest_stage}, elapsed {time.time() - t_cand:.1f}s)",
            flush=True,
        )

    rows.append(
        {
            "_meta": True,
            "kind": "headline",
            "candidates": summary,
            "writeback_to_catalogue": False,
            "elapsed_s": time.time() - t0,
        }
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")

    print(f"\n[#566] wrote {OUT_JSONL}", flush=True)
    print(f"[#566] total elapsed {time.time() - t0:.1f}s", flush=True)
    for s in summary:
        print(f"[#566] SUMMARY {s['candidate_id']}: {s['chain_verdict']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
