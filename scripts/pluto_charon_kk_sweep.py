"""Task #504: Pluto-Charon (k1,k2)-cycler family sweep — runnable script.

Sweeps (3,2) [positive control], (1,1), (2,1), (3,1), (2,2), (3,3) at
Pluto-Charon mu = cr3bp_system("Pluto","Charon").mu = 0.10876473603280369.
For each family: seeds from the nearest Ross-RT 2026 Table-I anchor, mu-steps
to PC mu, C-sweeps the family branch, finds the |nu|<1 stable window, and
brentq's the nu=0 midpoint.

Usage
-----
  uv run python scripts/pluto_charon_kk_sweep.py

Output: a per-family result table printed to stdout and saved to
  docs/notes/scratch/504_kk_sweep_raw.txt
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

# Make the src tree importable without `uv run` first installing the package.
_SRC = Path(__file__).parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from cyclerfinder.search.pluto_charon_kk_sweep import (  # noqa: E402
    PC_MU,
    SweepResult,
    make_pluto_charon_system,
    sweep_11,
    sweep_21,
    sweep_22,
    sweep_31,
    sweep_32_positive_control,
    sweep_33,
)


def _fmt(r: SweepResult) -> str:
    if r.stable_found:
        return (
            f"STABLE  C={r.jacobi_mid:.7f}  x0={r.x0_mid:.9f}  "
            f"T={r.period_mid:.5f} TU ({r.period_days:.3f} d)  "
            f"nu={r.nu_mid:.2e}  topo_ok={r.topology_ok}  xcheck={r.crosscheck_ok}"
        )
    return f"NEGATIVE  method={r.method!r}  note={r.note!r}"


def main() -> None:
    pc = make_pluto_charon_system()
    print(f"Pluto-Charon CR3BP sweep — mu={PC_MU}")
    print(f"  l_km={pc.l_km}  t_s={pc.t_s:.2f} s")
    print()

    sweeps = [
        ("(3,2) +ctrl", sweep_32_positive_control),
        ("(1,1)", sweep_11),
        ("(3,1)", sweep_31),
        ("(3,3)", sweep_33),
        ("(2,1)", sweep_21),
        ("(2,2)", sweep_22),
    ]

    lines: list[str] = []
    for label, fn in sweeps:
        print(f"  Running {label} ...", flush=True)
        t0 = time.time()
        result = fn()
        elapsed = time.time() - t0
        line = f"{label:14s}  {_fmt(result)}  [{elapsed:.1f}s]"
        print(f"  {line}")
        lines.append(line)

    print()
    print("Summary")
    print("-" * 80)
    for line in lines:
        print(line)

    out_path = Path(__file__).parent.parent / "docs" / "notes" / "scratch"
    out_path.mkdir(parents=True, exist_ok=True)
    (out_path / "504_kk_sweep_raw.txt").write_text("\n".join(lines) + "\n")
    print(f"\nSaved to {out_path / '504_kk_sweep_raw.txt'}")


if __name__ == "__main__":
    main()
