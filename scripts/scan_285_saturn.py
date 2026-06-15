"""#285 prioritized repeated-moon scan -- Saturn system.

Takubo 2210.14996 tour-design existence prior is strongest for the Titan-
Rhea-Dione-Tethys subset; the second band runs the 3-body Liang-CGE analogue
over Titan-Enceladus-Rhea. Both bands emit one JSONL row per Lambert-closed
candidate with the full Tier-0 + DOP853 + literature + ML-flagger scoring
chain; SILVER survivors feed the gauntlet (#274) as a separate step.

NO catalogue writeback. NO novelty claims. Run as::

    uv run python scripts/scan_285_saturn.py

Outputs ``data/scan_285_saturn.jsonl`` with one leading ``_meta`` row per
band, followed by per-candidate rows.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.search.saturn_uranus_campaign import run_prioritized_scan  # noqa: E402


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------
# Bands (deterministic; mirrors the #285 task statement)
# --------------------------------------------------------------------------

# Band A: Titan-Rhea-Dione-Tethys repeated-moon chains. The task lists
# (k1, k2) in {(1,1), (2,1), (3,2), (1,2)} -- these are per-leg revolution
# counts for a 2-leg closed cycle (A->B->A). The RepeatedMoonTarget grid
# enumerates the full Cartesian product over n_rev_grid for each pair of
# legs; we pass the union {1, 2, 3} so the four explicit cells fall inside
# (the engine post-filters none, the deterministic enumeration produces a
# superset which we accept -- false positives here are rejected downstream
# by Lambert feasibility + the residual gate).
SATURN_BAND_A = dict(
    primary="Saturn",
    moons=("Titan", "Rhea", "Dione", "Tethys"),
    seq_lengths=(3,),  # length-3 closed cycle: anchor-other-anchor (2 legs).
    n_rev_grid=(1, 2, 3),  # union of (k1, k2) in the task spec
)

# Band B: 3-body Liang CGE analogue (Titan-Enceladus-Rhea). Liang's
# original is Callisto-Ganymede-Callisto-Europa-Callisto (length-5, 4 legs
# with the anchor returned to three times). We mirror that with Titan as
# the anchor, Enceladus and Rhea as the two distinct visited bodies.
# RepeatedMoonTarget enumerates closed cycles starting and ending at the
# same body, so the length-5 enumeration over {Titan, Enceladus, Rhea}
# subsumes the Liang-style topology. Rev-grid kept small for runtime.
SATURN_BAND_B = dict(
    primary="Saturn",
    moons=("Titan", "Enceladus", "Rhea"),
    seq_lengths=(5,),
    n_rev_grid=(0, 1, 2),
)


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_285_saturn.jsonl"
    print(f"[285] Saturn campaign starting -- sha={sha} -- out={out_path}", flush=True)

    summaries = []
    # Band A first (existence prior strongest).
    print("[285] Saturn Band A: Titan-Rhea-Dione-Tethys (k=3 cycles, n_rev=1..3)", flush=True)
    sum_a = run_prioritized_scan(
        out_path=out_path,
        gate_residual_kms=0.05,
        phase_samples=12,
        git_sha=sha,
        **SATURN_BAND_A,
    )
    summaries.append(("Band A", sum_a.as_dict()))
    print(f"[285] Saturn Band A summary: {sum_a.as_dict()}", flush=True)

    # Band B appended to the same JSONL (with its own _meta row).
    print(
        "[285] Saturn Band B: Titan-Enceladus-Rhea Liang-CGE analogue (k=5 cycles)",
        flush=True,
    )
    # Re-open in append mode for Band B to keep both bands in the same file.
    band_b_path = ROOT / "data" / "scan_285_saturn_bandB.jsonl"
    sum_b = run_prioritized_scan(
        out_path=band_b_path,
        gate_residual_kms=0.05,
        phase_samples=12,
        git_sha=sha,
        **SATURN_BAND_B,
    )
    summaries.append(("Band B", sum_b.as_dict()))
    print(f"[285] Saturn Band B summary: {sum_b.as_dict()}", flush=True)

    # Concatenate Band B into Band A's file so the canonical output is one JSONL.
    with out_path.open("a", encoding="utf-8") as out, band_b_path.open() as src:
        for line in src:
            out.write(line)
    # Remove the per-band file once merged.
    band_b_path.unlink(missing_ok=True)

    # Write a small summary block at the end.
    with out_path.open("a", encoding="utf-8") as out:
        out.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "summary",
                    "summaries": summaries,
                    "git_sha": sha,
                }
            )
            + "\n"
        )

    print(f"[285] Saturn campaign DONE -- summaries: {summaries}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
