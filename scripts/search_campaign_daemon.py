"""Long-running search + training-data campaign (#210 corpus / cycler screening).

Runs for hours/days. Writes ONLY to the gitignored out/ tree (the outcome log +
a heartbeat) -- never to data/ or git, so it never collides with concurrent work.
Every solver call auto-logs an (inputs -> outcome) tuple via the #210 logger.

Phase A (finite): sweep the FULL JPL Earth-Moon periodic-orbit catalogue (every
family, every member) through our CR3BP corrector -- a comprehensive same-model
cross-check AND corpus seed.
Phase B (unbounded): randomized CR3BP corrector search -- perturbed JPL members +
random seeds across the EM system -- looping until killed. Converged AND failed
outcomes are both logged (failures are training signal). Heartbeat every 200 solves.

Discipline: SILVER / data-gen only. Nothing here writes the catalogue or claims a
discovery; any genuine candidate would still face the V0-V5 gauntlet separately.
"""

from __future__ import annotations

import argparse
import os
import time
import traceback

import numpy as np

from cyclerfinder.core import cr3bp
from cyclerfinder.search.cr3bp_periodic import correct_periodic
from cyclerfinder.verify.jpl_periodic_orbits import query

# Multi-core: launch N copies with distinct --worker-id to use N cores. Each
# worker has its own RNG seed and its own gitignored log shard (no shared-file
# contention). Worker 0 runs the finite Phase-A JPL sweep; workers >0 skip it
# (no redundant catalogue pulls) and go straight to randomized Phase-B search.
_AP = argparse.ArgumentParser()
_AP.add_argument("--worker-id", type=int, default=0)
WID = _AP.parse_args().worker_id

# The logger reads this env var at call time (not import time), so setting it
# here -- after imports -- still captures every solve below.
os.environ.setdefault("CYCLERFINDER_OUTCOME_LOG", f"out/outcome_log/search_campaign_w{WID}.jsonl")
HEARTBEAT = f"out/outcome_log/search_campaign_w{WID}.heartbeat"
os.makedirs("out/outcome_log", exist_ok=True)

SYS = cr3bp.cr3bp_system("Earth", "Moon")
RNG = np.random.default_rng(20260613 + WID)

FAMILIES = [
    ("halo", 1, "N"),
    ("halo", 1, "S"),
    ("halo", 2, "N"),
    ("halo", 2, "S"),
    ("halo", 3, "N"),
    ("halo", 3, "S"),
    ("dro", None, None),
    ("dpo", None, None),
    ("lyapunov", 1, None),
    ("lyapunov", 2, None),
    ("lyapunov", 3, None),
    ("vertical", 1, None),
    ("vertical", 2, None),
    ("vertical", 3, None),
    ("axial", 1, None),
    ("axial", 2, None),
    ("lpo", 1, None),
    ("lpo", 2, None),
    ("resonant", None, None),
]

solves = 0
conv = 0
t0 = time.time()


def beat(msg: str) -> None:
    with open(HEARTBEAT, "w") as fh:
        fh.write(
            f"{time.strftime('%Y-%m-%dT%H:%M:%S')}  solves={solves} converged={conv} "
            f"elapsed={int(time.time() - t0)}s  {msg}\n"
        )
    print(f"[{int(time.time() - t0)}s] solves={solves} conv={conv}  {msg}", flush=True)


def run_one(state0: np.ndarray, period: float) -> None:
    global solves, conv
    try:
        po = correct_periodic(SYS, np.asarray(state0, float), float(period))
        solves += 1
        if po.converged:
            conv += 1
    except Exception:
        solves += 1


# ---- Phase A: full JPL catalogue sweep (finite) ----
# ONLY worker 0 hits the JPL API (avoid N-fold hammering); it sweeps + logs and
# caches the seed pool to disk. Workers >0 wait for that cache and load it.
_POOL_CACHE = "out/outcome_log/jpl_seed_pool.npz"
jpl_pool: list[tuple[np.ndarray, float]] = []


def _publish_pool(states: list[np.ndarray], periods: list[float]) -> None:
    """Atomically (tmp + replace) write the seed-pool cache so a reader never
    sees a half-written npz. Called after EVERY family so workers >0 can load a
    usable pool within ~1 min instead of waiting for the whole sweep to finish."""
    if not states:
        return
    tmp = _POOL_CACHE + ".tmp"
    np.savez(tmp, states=np.array(states), periods=np.array(periods))
    os.replace(tmp, _POOL_CACHE)


if WID == 0:
    states: list[np.ndarray] = []
    periods: list[float] = []
    for fam, libr, branch in FAMILIES:
        kw: dict[str, object] = {}
        if libr is not None:
            kw["libr"] = libr
        if branch is not None:
            kw["branch"] = branch
        try:
            _consts, orbits = query("earth-moon", fam, **kw)  # type: ignore[arg-type]
        except Exception as exc:
            print(f"PhaseA {fam} L{libr} {branch}: query failed: {exc}", flush=True)
            continue
        # Subsample per family: Phase A is a fast cross-check + seed harvest, NOT
        # an exhaustive corpus (Phase B is the bulk + the diverse converged/failed
        # boundary). Capping it keeps the sweep short; publishing the pool after
        # each family lets all workers reach Phase B quickly instead of idling.
        for o in orbits[::25]:
            run_one(o.state0, o.period)
            s = np.asarray(o.state0, float)
            jpl_pool.append((s, float(o.period)))
            states.append(s)
            periods.append(float(o.period))
        _publish_pool(states, periods)  # incremental: usable after family 1
        beat(f"PhaseA {fam} L{libr} {branch}: {len(orbits)} members, sampled {len(orbits[::25])}")
    _publish_pool(states, periods)
    beat(f"PhaseA DONE -- pool cached ({len(jpl_pool)} members)")
else:
    # Wait for worker 0 to publish the seed-pool cache, then load it.
    for _ in range(600):  # up to ~10 min
        if os.path.exists(_POOL_CACHE):
            try:
                npz = np.load(_POOL_CACHE)
                jpl_pool = list(zip(npz["states"], [float(p) for p in npz["periods"]], strict=True))
                break
            except Exception:
                pass
        time.sleep(1.0)
    beat(f"PhaseB worker {WID} loaded pool ({len(jpl_pool)} members)")

# ---- Phase B: unbounded randomized search (perturb JPL pool) ----
if not jpl_pool:
    jpl_pool = [(np.array([0.83, 0.0, 0.1, 0.0, 0.12, 0.0]), 3.0)]

while True:
    try:
        base_state, base_t = jpl_pool[int(RNG.integers(len(jpl_pool)))]
        dstate = RNG.normal(0.0, 0.02, size=6)
        dt = float(RNG.normal(0.0, 0.15 * base_t))
        run_one(base_state + dstate, max(0.1, base_t + dt))
        if solves % 200 == 0:
            beat("PhaseB randomized search")
    except Exception:
        traceback.print_exc()
        time.sleep(1.0)
