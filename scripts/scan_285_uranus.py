"""#285 prioritized repeated-moon scan -- Uranus regulars (speculative scout).

No strong existence prior; Tisserand self-pruning is expected to leave the
band empty. The scout's value is the negative: a clean empty across all
three pairs + the 3-body case answers "are Uranus regulars within reach of
the repeated-moon multi-rev genome at all?" with the same scorer chain the
Saturn band uses.

Bands:

* Titania-Oberon, Titania-Umbriel, Oberon-Umbriel (2-body pairs, length-3 cycles)
* Ariel-Umbriel-Titania (3-body, length-5 Liang-CGE-style cycle)

NO catalogue writeback. NO novelty claims. Run as::

    uv run python scripts/scan_285_uranus.py

Outputs ``data/scan_285_uranus.jsonl``.
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


# Band A: 2-body Titania-Oberon, Titania-Umbriel, Oberon-Umbriel.
# RepeatedMoonTarget enumerates ALL closed length-3 cycles over the moon
# set (filtering open or degenerate ones), so passing the full
# {Titania, Oberon, Umbriel} triple subsumes the three explicit pairs;
# the campaign's deterministic enumeration walks every (A, B, A) and
# (B, A, B) combination across the triple.
URANUS_BAND_A = dict(
    primary="Uranus",
    moons=("Titania", "Oberon", "Umbriel"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3),
)

# Band B: 3-body Ariel-Umbriel-Titania length-5 cycles.
URANUS_BAND_B = dict(
    primary="Uranus",
    moons=("Ariel", "Umbriel", "Titania"),
    seq_lengths=(5,),
    n_rev_grid=(0, 1, 2),
)


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_285_uranus.jsonl"
    print(f"[285] Uranus campaign starting -- sha={sha} -- out={out_path}", flush=True)

    summaries = []
    print("[285] Uranus Band A: Titania-Oberon-Umbriel pairs (k=3, n_rev=0..3)", flush=True)
    sum_a = run_prioritized_scan(
        out_path=out_path,
        gate_residual_kms=0.05,
        phase_samples=12,
        git_sha=sha,
        **URANUS_BAND_A,
    )
    summaries.append(("Band A", sum_a.as_dict()))
    print(f"[285] Uranus Band A summary: {sum_a.as_dict()}", flush=True)

    print("[285] Uranus Band B: Ariel-Umbriel-Titania Liang-CGE analogue (k=5)", flush=True)
    band_b_path = ROOT / "data" / "scan_285_uranus_bandB.jsonl"
    sum_b = run_prioritized_scan(
        out_path=band_b_path,
        gate_residual_kms=0.05,
        phase_samples=12,
        git_sha=sha,
        **URANUS_BAND_B,
    )
    summaries.append(("Band B", sum_b.as_dict()))
    print(f"[285] Uranus Band B summary: {sum_b.as_dict()}", flush=True)

    with out_path.open("a", encoding="utf-8") as out, band_b_path.open() as src:
        for line in src:
            out.write(line)
    band_b_path.unlink(missing_ok=True)

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

    print(f"[285] Uranus campaign DONE -- summaries: {summaries}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
