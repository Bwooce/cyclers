"""Cross-check our CR3BP corrector against the JPL Three-Body Periodic Orbit
Catalog (#116 oracle). NETWORK script (live JPL query) — not a test.

For a sample of JPL-published Earth-Moon halo ICs we:
  1. reconcile the mass-ratio convention (JPL vs ours);
  2. re-propagate each JPL IC for the JPL-published period under OUR integrator,
     under BOTH mu values, and report the periodicity-closure residual
     |X(T)-X(0)| (the JPL IC is a residual-zero point in their model; how well
     it re-closes in ours, and how much of any gap is the mu mismatch, is the
     per-interface external anchor the false-consensus discipline wants);
  3. compare OUR Jacobi constant at the IC against JPL's published Jacobi.

This certifies our CR3BP propagator + Jacobi against an independent group's
catalog. It is a verification artifact, NOT a catalogue writeback (halo orbits
are not cyclers). Run: ``uv run python scripts/jpl_oracle_crosscheck.py``.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core import cr3bp
from cyclerfinder.verify.jpl_periodic_orbits import query, reconcile_mu


def _closure_residual(system: cr3bp.CR3BPSystem, state0: np.ndarray, period: float) -> float:
    arc = cr3bp.propagate(system, state0, period, with_stm=False)
    return float(np.linalg.norm(arc.state_f - state0))


def main() -> None:
    print("JPL periodic-orbit oracle cross-check (Earth-Moon halo)")
    print("=" * 64)
    constants, orbits = query("earth-moon", "halo", libr=2, branch="N")
    rec = reconcile_mu(constants.mu)
    print(f"JPL mu = {rec['jpl_mu']:.15g}")
    print(f"our mu = {rec['our_mu']:.15g}")
    print(f"abs diff = {rec['abs_diff']:.3e}   rel diff = {rec['rel_diff']:.3e}")
    print(f"JPL L2-halo-N orbits returned: {len(orbits)}")
    print()

    # Our system at our mu, and a JPL-mu twin, so we can separate "our integrator
    # disagrees" from "the mu convention differs".
    sys_ours = cr3bp.cr3bp_system("Earth", "Moon")
    # Same dynamics, JPL's mu (length/time units irrelevant to the nd closure).
    sys_jpl = cr3bp.CR3BPSystem(
        mu=constants.mu,
        primary="Earth",
        secondary="Moon",
        l_km=sys_ours.l_km,
        t_s=sys_ours.t_s,
    )

    # Sample across the family (first, quartiles, last) rather than just the top.
    n = len(orbits)
    idxs = sorted({0, n // 4, n // 2, 3 * n // 4, n - 1})
    print(
        f"{'idx':>5} {'JPL C':>12} {'our C@ourmu':>14} {'dC':>10} "
        f"{'close@ourmu':>12} {'close@jplmu':>12}"
    )
    print("-" * 70)
    worst_jplmu = 0.0
    for i in idxs:
        o = orbits[i]
        our_c = cr3bp.jacobi_constant(o.state0, sys_ours.mu)
        dc = our_c - o.jacobi
        r_ours = _closure_residual(sys_ours, o.state0, o.period)
        r_jpl = _closure_residual(sys_jpl, o.state0, o.period)
        worst_jplmu = max(worst_jplmu, r_jpl)
        print(f"{i:>5} {o.jacobi:>12.6f} {our_c:>14.6f} {dc:>10.2e} {r_ours:>12.2e} {r_jpl:>12.2e}")
    print()
    print(f"worst closure residual under JPL mu: {worst_jplmu:.3e} (nd)")
    print("Reading: residual@jplmu isolates integrator agreement (mu matched);")
    print("the ours-vs-jpl gap at our mu is dominated by the ~1e-7 mu offset.")


if __name__ == "__main__":
    main()
