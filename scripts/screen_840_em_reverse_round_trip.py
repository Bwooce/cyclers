"""#840: Vaquero 3:1 -> 2:1 Earth-Moon reverse-direction free transfer at the
two CATALOGUED band-edge Jacobi constants (C=2.54, C=2.66).

`#822` demonstrated the reverse direction Wu(3:1)->Ws(2:1) at C=2.60 only -- a
C at which neither family node is a catalogued row. `#828` recorded the
FORWARD direction Wu(2:1)->Ws(3:1) touching two catalogued rows at their own
band edges (`vaquero-31-c254-em-cycler-2013` at C=2.54,
`vaquero-21-c266-em-cycler-2013` at C=2.66). This script runs the REVERSE
direction at those SAME two C values, at this project's own registry mu
(`cyclerfinder.core.cr3bp.cr3bp_system("Earth", "Moon").mu`), so each row ends
up with a demonstrated round trip (both directions, same C, same mu).

METHOD (this task's own diagnosed result, not the naive approach): a first
attempt via the plain blind
:func:`~cyclerfinder.search.vaquero_em_cycler_connections.find_free_transfer_reverse`
scan (the direct `#822` Sec. 4 template) at the sibling-module default
``epsilon=1e-4`` converged 6/11 refined candidates at C=2.66, but EVERY ONE
failed ``verify_connection``'s forward-reapproach gate -- including the seed
closest to the CR3BP time-reversal-symmetry-predicted crossing itself
(``forward_distance=1.58`` vs the ``0.5`` ceiling; module docstring, "the
same evidence-quality knob `#822`'s own docstring documents"). This script
instead uses
:func:`~cyclerfinder.search.vaquero_em_cycler_connections.find_free_transfer_reverse_near_prediction`,
which (a) predicts each C's reverse crossing directly from the ALREADY-
VERIFIED forward connection's own crossing via
:func:`~cyclerfinder.search.vaquero_em_cycler_connections.predict_reverse_crossing`
(no search needed for the prediction itself -- both orbits are x-axis
symmetric), and (b) sweeps epsilon upward, trying up to
``max_try_per_epsilon`` predicted-ranked candidates per epsilon, when the
closest-to-predicted candidate converges but fails verify_connection's gates
(never loosening the gates themselves). Result: C=2.54 passes immediately at
epsilon=1e-4 (unchanged) on the very first/closest-to-predicted candidate --
its converged crossing matches the symmetry prediction to ~6e-9, confirming
it is genuinely the time-reversal image of the SAME forward connection. C=2.66
needs epsilon=2e-4, and the candidate that verifies there is a DIFFERENT
intersection point (~0.26 away from the symmetry prediction, not its mirror
image) -- still a genuine, independently verified Wu(3:1)->Ws(2:1) connection
at the SAME C and mu (`#822`'s own module docstring already documents this
family pair's manifold intersection structure as "rich", not a single point).

The forward crossings are read directly from
``data/found/822_vaquero_em_free_transfer/results.json`` (never
hand-transcribed) so this script cannot silently drift from that record.

Writes the full per-C record (connection parameters + complete verification
evidence) to ``data/found/840_em_reverse_round_trip/results.json``.

Foreground, single-process. Observed runtime this task: ~1 min at C=2.54
(passes on the first epsilon), ~10 min at C=2.66 (two epsilon rounds, several
candidates each).
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.vaquero_em_cycler_connections as vcc

_REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = _REPO_ROOT / "data" / "found" / "840_em_reverse_round_trip"
FORWARD_RESULTS_PATH = (
    _REPO_ROOT / "data" / "found" / "822_vaquero_em_free_transfer" / "results.json"
)

#: The two catalogued band-edge Jacobi constants (#828's two row-touching entries).
C_POINTS = [2.54, 2.66]


def _forward_crossings() -> dict[float, tuple[float, float]]:
    """Each C's forward Wu(2:1)->Ws(3:1) crossing, read from #822's own
    sweep record (never hand-transcribed)."""
    raw = json.loads(FORWARD_RESULTS_PATH.read_text(encoding="utf-8"))
    out: dict[float, tuple[float, float]] = {}
    for row in raw["sweep"]:
        c = round(row["jacobi"], 4)
        if c in C_POINTS and row["verified"]:
            xv = row["connection"]["crossing_xv"]
            out[c] = (float(xv[0]), float(xv[1]))
    return out


def _conn_dict(conn: Any) -> dict[str, Any]:
    d = asdict(conn)
    d["crossing_xv"] = [float(v) for v in conn.crossing_xv]
    return d


def _ev_dict(ev: Any) -> dict[str, Any]:
    d = asdict(ev)
    d["crossing_state"] = [float(v) for v in ev.crossing_state]
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--c-values",
        nargs="*",
        type=float,
        default=None,
        help="subset of C points to (re)run this invocation (default: C_POINTS)",
    )
    args = ap.parse_args()

    system = cr3bp.cr3bp_system("Earth", "Moon")
    assert abs(system.mu - 0.01215058439469525) < 1e-15, f"registry mu drifted: {system.mu!r}"

    forward = _forward_crossings()
    todo = args.c_values if args.c_values is not None else C_POINTS

    rows: list[dict[str, Any]] = []
    for c in todo:
        c = round(float(c), 4)
        predicted = vcc.predict_reverse_crossing(forward[c])
        print(
            f"[{datetime.now(UTC).isoformat(timespec='seconds')}] C={c}: predicted reverse "
            f"crossing (x, xdot) = {predicted} ...",
            flush=True,
        )
        res = vcc.find_free_transfer_reverse_near_prediction(system, c, predicted)
        row: dict[str, Any] = {
            "jacobi": c,
            "direction": "Wu(3:1)->Ws(2:1)",
            "method": "find_free_transfer_reverse_near_prediction",
            "predicted_crossing_xv": list(predicted),
            "mu": system.mu,
            "n_seeds": res.n_seeds,
            "n_refined": res.n_refined,
            "n_converged": res.n_converged,
            "verified": res.connection is not None,
            "notes": res.notes,
        }
        if res.connection is not None and res.evidence is not None:
            conn_d = _conn_dict(res.connection)
            ev_d = _ev_dict(res.evidence)
            row["connection"] = conn_d
            row["evidence"] = ev_d
            print(
                f"  C={c}: VERIFIED epsilon={conn_d['epsilon']:g} "
                f"k=({conn_d['k_u']},{conn_d['k_s']}) "
                f"branches=({conn_d['branch_u']:+d},{conn_d['branch_s']:+d}) "
                f"crossing=({conn_d['crossing_xv'][0]:+.10f},{conn_d['crossing_xv'][1]:+.10f}) "
                f"residual={conn_d['residual']:.3e} full_gap={ev_d['full_state_gap']:.3e} "
                f"radau_gap={ev_d['radau_gap']:.3e} "
                f"back={ev_d['backward_distance']:.3e} fwd={ev_d['forward_distance']:.3e}",
                flush=True,
            )
        else:
            print(f"  C={c}: NOT VERIFIED -- {res.notes}", flush=True)
        rows.append(row)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "results.json"
    out = {
        "task": "#840",
        "generated_utc": datetime.now(UTC).isoformat(),
        "mu": system.mu,
        "l_km": system.l_km,
        "t_s": system.t_s,
        "method_note": (
            "Symmetry-informed search (find_free_transfer_reverse_near_prediction): "
            "the blind find_free_transfer_reverse scan at epsilon=1e-4 (sibling-module "
            "default) converged 6/11 refined candidates at C=2.66 but every one failed "
            "verify_connection's forward-reapproach gate, including the "
            "closest-to-predicted candidate itself (forward_distance=1.58 vs the 0.5 "
            "ceiling). Raising epsilon (the same documented #822 evidence-quality knob) "
            "and trying multiple predicted-ranked candidates per epsilon "
            "(max_try_per_epsilon=5) found a verified connection at epsilon=2e-4 for "
            "C=2.66 (1e-4 sufficed unchanged at C=2.54, on the closest-to-predicted "
            "candidate). C=2.54's verified crossing matches the symmetry prediction to "
            "~6e-9 (genuinely the time-reversal image of the forward connection); "
            "C=2.66's verified crossing is a DIFFERENT intersection point (~0.26 away "
            "from the prediction) -- still a genuine, independently verified connection "
            "at the same C and mu. Gates were never loosened."
        ),
        "sweep": rows,
    }
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")

    n_ok = sum(1 for r in rows if r["verified"])
    print(f"\nverified {n_ok}/{len(rows)} C points; results at {path}", flush=True)


if __name__ == "__main__":
    main()
