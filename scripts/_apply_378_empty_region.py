"""Apply the #378 cislunar-BCT clean-negative empty-region record (Phase 4.1).

One-shot writer (prefixed `_` so it is treated as scratch and left untouched by
the harness): appends the WSB/BCT clean-negative EmptyRegionReport to
data/empty_regions.jsonl with a method-capability re-open key so a future
coherent-QBCP or stronger method subsumes it.
"""

from __future__ import annotations

import subprocess

from cyclerfinder.data.empty_regions import (
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability

git_sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()

method = MethodCapability(
    genome="#378 cislunar BCT substrate -- Belbruno WSB surface (core/wsb.py: E_2, "
    "periapsis, analytic-W eq 3.29, numerical stability-class) + backward BCT "
    "constructor (genome/bct_transfer.construct_bct_backward) on the incoherent "
    "BCR4BP (core/bcr4bp, Andreu/Rosales-Jorba constants). theta_2 family sweep.",
    corrector="forward 2x2 (|V0|,gamma0)->(r_23, E_2-on-W) terminal targeting "
    "(genome/bct_transfer.correct_bct_forward); on-W ballistic capture from LEO "
    "did NOT converge (documented incoherent-model gap, design R3).",
    capability_tags=frozenset(
        {
            "bcr4bp",
            "incoherent-bcr4bp",
            "wsb-surface",
            "ballistic-capture",
            "bct-backward-construction",
            "cislunar",
            "sun-earth-moon",
            "theta2-family-sweep",
            "return-leg-w-reacquisition",
        }
    ),
    git_sha=git_sha,
)

report = EmptyRegionReport(
    region_id="cislunar-bct-wsb-quasicycler-2026-06-26",
    family="Sun-Earth-Moon cislunar weak-stability-boundary (WSB) ballistic-capture "
    "transfers (Belbruno 2004 Ch.3) -- searched for a WSB-anchored quasi-cycler: a "
    "repeating capture<->escape chain whose return leg re-acquires W (the only "
    "catalogue-class candidate; a single BCT is a precursor_mga transfer, not a cycler).",
    centre="Earth-Moon (mu=0.01215) BCR4BP with the Andreu Sun perturbation; "
    "Hiten-class exterior WSB-transfer regime (Sun-shaped Earth apoapsis ~3.9 LD).",
    topologies=(
        {
            "object": "exterior WSB transfer (Hiten-class)",
            "capture": "ballistic (E_2<=0) at r_M+100km, e_2~0.95",
            "apoapsis_band_ld": [2.7, 5.1],
            "theta2_family": "0.50..1.60 rad",
        },
    ),
    method_capability=method,
    search_extent={
        "points_total": 44,
        "theta2_values": [0.50, 0.60, 0.65, 0.70, 0.75, 0.80, 1.00, 1.05, 1.25, 1.40, 1.60],
        "e2_values": [0.90, 0.95],
        "branches": ["retrograde", "direct"],
        "back_days": 70.0,
        "return_leg_forward_days": 90.0,
    },
    prune_gates=(
        "max_earth_apoapsis_ld in [2.7,5.1] (Hiten signature band)",
        "capture E_2<=0 (on W, ballistic)",
        "return_leg_reacquires_W (E_2<=0 AND periapsis AND r<5*r0 after leaving capture)",
    ),
    result={
        "n_evaluated": 44,
        "n_in_hiten_apoapsis_band": 4,
        "best_apoapsis_ld": 3.95,
        "best_apoapsis_theta2": 1.25,
        "n_quasi_cycler_candidates": 0,
        "capture_e2_on_w": True,
        "forward_from_leo_on_w_capture_converged": False,
    },
    verdict="EMPTY",
    interpretation="CAPABILITY-ONLY, no catalogue-class object. The WSB surface + "
    "backward BCT constructor reproduce the Hiten apoapsis signature (3.95 LD, "
    "bullseye on Belbruno's 3.9 LD) with exact on-W ballistic capture (E_2<0), but "
    "NO BCT's return leg re-acquires W -- no WSB-anchored cislunar quasi-cycler. "
    "Consistent with Belbruno Theorem 3.58 (capture on W is a CHAOTIC process, "
    "cutting against clean periodicity). Two documented model-gap boundaries remain "
    "re-open keys: (1) the forward-from-LEO on-W capture did not converge in the "
    "incoherent BCR4BP, and (2) the coherent QBCP (alpha-table, un-digested) is a "
    "higher-fidelity model. A repeating object, if any exists, needs one of these.",
    source_anchors="Belbruno 2004 (Capture Dynamics; Def 3.11 ballistic capture, "
    "Def 3.12/eq 3.9 WSB, Lemma 3.21/eq 3.29 analytic W, §3.4 Hiten BCT, Thm 3.58 "
    "chaotic capture); Koon-Lo-Marsden-Ross 2001 (Shoot the Moon manifold mechanism). "
    "digest: docs/notes/2026-06-17-digest-belbruno-2004.md.",
    run={
        "date": "2026-06-26",
        "git_sha": git_sha,
        "issue": 378,
        "grid_points": 44,
        "note": "docs/notes/2026-06-26-378-cislunar-bct-verdict.md",
    },
)

append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
print(f"appended cislunar-BCT clean-negative to {DEFAULT_EMPTY_REGIONS_PATH} (sha {git_sha})")
