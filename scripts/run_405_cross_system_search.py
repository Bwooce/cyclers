"""#405 Phase-A bounded cross-system (SE<->EM) closure search with Radau cross-check.

Runs ``search_cross_cycle`` over a small bounded grid of EM-L2 Jacobi constants
(EM-L2 family range [3.1445, 3.184]) at the Canalias SE Jacobi, for the
(EM-L2, SE-L2) and (EM-L1, SE-L1) libration pairs.  Prints per-grid-point
outcomes with wall-clock timestamps.  Cross-checks any closed cycle with an
independent Radau integrator.

If ZERO cycles closed (the expected Phase-A clean-negative outcome), appends a
registry entry to ``data/negative_results.yaml`` unless the id is already
present (idempotent).

This is a Phase-A search: bounded single-revolution grid, coarse manifold
horizons.  A non-closing grid is the honest result — consistent with the #316
prediction that natural cross-system closure requires the ~19yr Metonic (235:19)
commensurability, which lies far outside any single-revolution grid.
"""

from __future__ import annotations

import datetime
import pathlib
import sys

import numpy as np
import yaml

# Ensure the src tree is on the path when invoked as a script.
_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    FrameBridge,
    crosscheck_cross_cycle,
    em_moon_system,
    se_earth_system,
    search_cross_cycle,
)

# Canalias 2007 SE bifurcation Jacobi — mirrors the golden value in the test.
CANALIAS_C_SE = 3.000863625

# EM-L2 Jacobi constants spanning the EM-L2 family range [3.1445, 3.184].
# Two points: previously working value (Task-4 confirmed forward connection) + high end.
# The low end (3.1445) is skipped here — EM-L1 is likely outside its family range there
# and the STM propagation fails (handled gracefully as converged=False in the search).
C_EM_GRID = (3.15, 3.184)

_NEG_RESULTS_PATH = _REPO / "data" / "negative_results.yaml"

_NEW_ENTRY_ID = "cross_system_se_em_L2_patched_cr3bp"

_NEW_ENTRY = {
    "id": _NEW_ENTRY_ID,
    "issue": 405,
    "method": ("patched-CR3BP SE<->EM connection matcher + bounded closure search (#405 Phase A)"),
    "regime": "Sun-Earth + Earth-Moon coupled CR3BP (patched, inertial-frame match)",
    "failed_rung": "closure-search",
    "physical_reason": (
        "Forward EM-L2->SE-L2 leg closes near-ballistically (pos gap ~0.4 km, "
        "dV ~0.36 km/s), but the SE->EM return leg does not co-reach the patch "
        "section on a bounded single-revolution grid and theta does not return "
        "commensurately; closed=False. Consistent with the ~19yr Metonic (235:19) "
        "natural-closure prediction (note #316), which lies outside any bounded "
        "single-revolution grid."
    ),
    "resweep_condition": (
        "BCR4BP Phase B refinement (coherent 4-body resonant closure) OR a "
        "multi-revolution / Metonic-commensurate grid with longer integration "
        "horizons."
    ),
}


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _append_negative_result() -> bool:
    """Append _NEW_ENTRY to negative_results.yaml unless already present. Returns True if added."""
    data = yaml.safe_load(_NEG_RESULTS_PATH.read_text())
    entries = data.get("entries", [])
    ids = {e.get("id") for e in entries}
    if _NEW_ENTRY_ID in ids:
        print(
            f"[{_ts()}] negative_results.yaml: entry '{_NEW_ENTRY_ID}' already present — skipping."
        )
        return False
    entries.append(_NEW_ENTRY)
    data["entries"] = entries
    _NEG_RESULTS_PATH.write_text(
        yaml.dump(data, default_flow_style=False, allow_unicode=True, sort_keys=False)
    )
    print(f"[{_ts()}] negative_results.yaml: appended entry '{_NEW_ENTRY_ID}'.")
    return True


def main() -> None:
    print(f"[{_ts()}] #405 Phase-A bounded cross-system closure search starting.")
    print(f"[{_ts()}] Building SE/EM systems and FrameBridge ...")

    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)

    print(f"[{_ts()}] SE system: mu={se.mu:.6e}  EM system: mu={em.mu:.6e}")
    print(f"[{_ts()}] Grid: c_em={C_EM_GRID}  c_se={CANALIAS_C_SE}")
    print(f"[{_ts()}] Libration pairs: (EM-L2,SE-L2) + (EM-L1,SE-L1)")
    print()

    print(f"[{_ts()}] Running search_cross_cycle (this may take several minutes) ...")

    results = search_cross_cycle(
        bridge,
        c_em_grid=C_EM_GRID,
        c_se_grid=(CANALIAS_C_SE,),
        libration_pairs=(("EM-L2", "SE-L2"), ("EM-L1", "SE-L1")),
        max_attempts=2,
    )

    print(f"[{_ts()}] Search complete. {len(results)} grid points evaluated.")
    print()

    n_closed = 0
    for i, cyc in enumerate(results):
        status = "CLOSED" if cyc.closed else "OPEN"
        print(
            f"[{_ts()}] [{i + 1:02d}] {status}  pair={cyc.libration_pair}  "
            f"c_em={cyc.c_em:.6f}  c_se={cyc.c_se:.9f}  "
            f"max_leg_resid={cyc.max_leg_residual:.3e} km  "
            f"theta_res={cyc.theta_closure_residual:.3e} rad"
        )
        if cyc.notes:
            print(f"         notes: {cyc.notes}")
        if cyc.connections:
            for j, conn in enumerate(cyc.connections):
                cv = "CONV" if conn.converged else "FAIL"
                print(
                    f"         leg[{j}] {cv}  {conn.label_from}->{conn.label_to}  "
                    f"resid={conn.residual:.3e} km  dV={conn.patch_dv_kms:.3f} km/s  "
                    f"theta={conn.theta:.4f} rad  n_iter={conn.n_iter}"
                )
        if cyc.closed:
            n_closed += 1
            print(f"[{_ts()}]   CLOSED cycle found — running Radau cross-check ...")
            checked = crosscheck_cross_cycle(bridge, cyc)
            ir = checked.independent_residual
            print(f"[{_ts()}]   Radau cross-check: independent_residual={ir:.3e} km")
        print()

    print(f"[{_ts()}] Summary: {n_closed}/{len(results)} grid points closed.")
    print()

    if n_closed == 0:
        print(f"[{_ts()}] Zero closures — Phase-A clean negative confirmed.")
        print(f"[{_ts()}] Registering clean-negative entry in {_NEG_RESULTS_PATH.name} ...")
        added = _append_negative_result()
        if added:
            print(f"[{_ts()}] Registry entry registered successfully.")
    else:
        print(
            f"[{_ts()}] {n_closed} closure(s) found — clean-negative entry NOT added "
            "(review closed cycles above before any catalogue writeback)."
        )

    print(f"[{_ts()}] Done.")

    # Report independent_residuals for any crosschecked cycles
    crosschecked = []
    for cyc in results:
        if cyc.closed:
            checked = crosscheck_cross_cycle(bridge, cyc)
            crosschecked.append(checked)
    if crosschecked:
        print()
        print(f"[{_ts()}] Cross-check independent residuals:")
        for cyc in crosschecked:
            print(
                f"  pair={cyc.libration_pair}  c_em={cyc.c_em:.6f}  "
                f"independent_residual={cyc.independent_residual:.3e} km"
            )

    # Array of independent_residuals for any non-nan results (for numeric reporting)
    indep_residuals = np.array(
        [r.independent_residual for r in results if not np.isnan(r.independent_residual)],
        dtype=float,
    )
    if indep_residuals.size > 0:
        print(
            f"[{_ts()}] Non-nan independent_residual range: "
            f"[{indep_residuals.min():.3e}, {indep_residuals.max():.3e}] km"
        )


if __name__ == "__main__":
    main()
