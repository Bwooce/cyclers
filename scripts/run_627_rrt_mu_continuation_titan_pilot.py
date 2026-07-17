"""Task #627: Ross-Roberts-Tsoukkas 2026 mu-continuation PILOT to Saturn-Titan.

#494 continued the Ross & Roberts-Tsoukkas 2026 (arXiv:2606.29189) stable
ballistic prograde (k1,k2)-cycler families UP in mass ratio mu, from Earth-Moon
(mu=0.01215) toward real binary systems (Pluto-Charon mu=0.1085 through
mu=0.5), using a solver #494 itself characterizes as mu-agnostic. Nobody has
continued DOWN below the paper's own mu=0.001 floor toward a REAL planet-moon
system, all of which sit well below that floor (Saturn-Titan mu=2.37e-4,
Neptune-Triton ~2.1e-4, Jupiter-Ganymede ~7.8e-5, Jupiter-Europa ~2.5e-5,
Uranus-Titania ~4e-5).

THIS IS A PILOT ONLY (task #627, dispatched from the #623 strategic review
shortlist B4) -- continue exactly two representative family members, k=(1,1)
and k=(3,3), down to ONE target (Saturn-Titan) and check two gates:

  (a) STABILITY: Barden |nu| < 1 at the landed Titan-mu member.
  (b) ENCOUNTER-RELEVANCE: is the perimoon passage a genuine flyby-scale
      encounter (periapsis altitude + relative speed within a sane range,
      objectively judged via the existing physical_sanity.py bend gate) or a
      degenerate near-collision / too-distant pass?

Honest framing (see the #627 OUTSTANDING.md bullet and
docs/notes/2026-07-17-623-strategic-review.md section B4): as mu -> 0,
symmetric stable prograde periodic orbits shade into classical Poincare
first-kind/resonant families -- existence at small mu is close to
theorem-guaranteed. Surviving down to Titan mu is EXPECTED and is NOT by
itself evidence of anything novel. A literature novelty check
(search/literature_check.py) is mandatory before anything here is called a
discovery (see the companion notes file for that check's live-WebSearch
results; this script also runs the deterministic offline-corpus backend for
the record).

Continuation-machinery finding (the actual engineering deliverable of this
pilot): #494/#549's ADMITTED real-binary extension used a fixed-Jacobi
natural-parameter mu-step (:func:`real_binary_kk_sweep.mu_step_to_system`,
40 linear steps, holding C fixed the whole way) because going UP in mu never
crosses below the anchor's own C_L1. Going DOWN toward Titan mu is different:
C_L1(Saturn-Titan mu=2.37e-4) = 3.0158, which is BELOW both the (1,1) and
(3,3) Table-I mu=0.01215 anchors' own Jacobi constants (3.1512 / 3.1834) --
so a fixed-C walk silently drifts onto a topologically wrong,
non-secondary-reaching branch well before reaching the target (verified: the
naive mu_step_to_system "converges" at every step but lands on a (k1,0)
orbit that never leaves the primary's realm). This pilot instead uses
:func:`real_binary_kk_sweep.mu_step_to_system_tracking_c_l1`, a new function
this task adds that interleaves small C-walks (tracking just below C_L1(mu))
with each mu step -- the direct generalization of #494's own "Strategy B"
fix (``pluto_charon_kk_sweep.sweep_31``) for the same problem going up.

No catalogue writeback (this is a pilot; see task doc step 8).

Usage
-----
  uv run python scripts/run_627_rrt_mu_continuation_titan_pilot.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

_SRC = Path(__file__).parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np  # noqa: E402

import cyclerfinder.search.cr3bp_periodic as cp  # noqa: E402
from cyclerfinder.core.cr3bp import cr3bp_system  # noqa: E402
from cyclerfinder.core.satellites import SATELLITES  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.search.binary_star_search import winding_topology  # noqa: E402
from cyclerfinder.search.literature_check import (  # noqa: E402
    KNOWN_CORPUS,
    CandidateSignature,
    SearchResult,
    check_literature,
)
from cyclerfinder.search.mu_continuation import scan_c_family_at_mu  # noqa: E402
from cyclerfinder.search.perimoon_passage import find_perimoon_passage  # noqa: E402
from cyclerfinder.search.physical_sanity import flyby_is_useful  # noqa: E402
from cyclerfinder.search.pluto_charon_kk_sweep import _c_l1  # noqa: E402
from cyclerfinder.search.real_binary_kk_sweep import (  # noqa: E402
    ANCHORS,
    mu_step_to_system_tracking_c_l1,
)

OUT_PATH = Path(__file__).parent.parent / "docs" / "notes" / "scratch" / "627_titan_pilot_raw.txt"

_REGION_ID = "rrt-mu-continuation-titan-pilot-2026-07-17"
_METHOD = MethodCapability(
    genome=(
        "Ross-Roberts-Tsoukkas 2026 (k1,k2)-cycler CR3BP genome (#494's fixed-Jacobi symmetric "
        "corrector + Barden stability + winding-topology classifier + independent Radau "
        "crosscheck), continued DOWN in mass ratio mu from the Earth-Moon Table-I anchors "
        "toward Saturn-Titan mu=2.37e-4 (below the paper's own mu=0.001 floor)"
    ),
    corrector=(
        "mu_step_to_system_tracking_c_l1 (#627): interleaved C-walk-below-C_L1(mu) + "
        "fixed-Jacobi mu-step, generalizing #494's sweep_31 Strategy B to every step of a "
        "downward mu walk"
    ),
    capability_tags=frozenset(
        {"cr3bp", "binary-cycler", "k1k2-genome", "real-moon", "mu-continuation-downward"}
    ),
    git_sha="working-tree",
)

# (1,1): the mu=0.001 Table-I anchor (Rep 1) is used as the starting point, not the
# mu=0.01215 Rep-2 anchor the task's own bullet names -- #627 found empirically that
# the Rep-2 (1,1) member sits at (or immediately adjacent to) a fold: even a dmu=1e-7
# step off that exact anchor jumps the Newton corrector onto a DIFFERENT, unrelated
# branch (x0 discontinuously jumps -0.768 -> -0.512). The mu=0.001 anchor is itself a
# sourced Table-I row (still traces to Ross-RT 2026, not a fresh guess) and is free of
# this fold, so it is used as an intermediate, already-validated waypoint -- consistent
# with "however many intermediate continuation steps are numerically prudent" per the
# task's own step-size-discipline instruction.
REPRESENTATIVES: dict[str, dict] = {
    "11": dict(
        label="(1,1)",
        anchor_mu=0.001,
        x0=-0.647047499999966,
        jacobi=3.031605708907296,
        period=14.774502790974823,
        hc=1,  # empirically determined (nearest y=0 crossing to T/2 at this anchor)
        k1=1,
        k2=1,
        rrt_anchor_note="Table I Rep 1 (mu=0.001) -- the paper's own floor",
        # Validated in this task's own exploration at c_margin=0.02 (a ~7-minute run,
        # reproduced verbatim by tests/search/test_627_titan_pilot.py's slow test);
        # kept identical here so the script's own run matches that already-checked
        # result exactly.
        c_margin=0.02,
    ),
    "33": dict(
        label="(3,3)",
        anchor_mu=ANCHORS["mu01215_33"]["anchor_mu"],
        x0=ANCHORS["mu01215_33"]["x0"],
        jacobi=ANCHORS["mu01215_33"]["jacobi"],
        period=ANCHORS["mu01215_33"]["period"],
        hc=7,  # empirically determined (nearest y=0 crossing to T/2 at this anchor)
        k1=3,
        k2=3,
        rrt_anchor_note="Table I Rep 3 (mu=0.01215) -- the task's literal starting point",
        # 0.02 is too large here: the (3,3) anchor's own C already sits only 0.0049
        # below C_L1(0.01215), so a margin of 0.02 forces an immediate C-walk even at
        # the FIRST (zero-mu-change) step, and the hc=7 branch is not robust to a
        # C-decrease that large at fixed mu (found empirically -- #627). 0.005 avoids
        # forcing that walk while still tracking C_L1(mu) down the rest of the range.
        c_margin=0.005,
    ),
}


def _append(lines: list[str]) -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("a") as f:
        f.write("\n".join(lines) + "\n")


def _print_and_append(lines: list[str]) -> None:
    for line in lines:
        print(line, flush=True)
    _append(lines)


def offline_corpus_search(query: str) -> list[SearchResult]:
    """Deterministic offline literature backend (mirrors #494/#582's own use).

    Only re-finds families already in the curated corpus; kept here for the
    record alongside the live-WebSearch check documented in this task's notes
    file (a structural-fingerprint offline backend cannot itself perform a
    real literature search -- see module docstring of literature_check.py).
    """
    q = query.lower()
    out: list[SearchResult] = []
    for anchor in KNOWN_CORPUS:
        hit = any(a.lower() in q for a in anchor.authors) or any(
            kw.lower() in q for kw in anchor.keywords
        )
        bodies_named = sum(1 for b in anchor.body_set if b.lower() in q)
        if hit or (bodies_named >= 2 and "cycler" in q):
            bodies = " ".join(sorted(anchor.body_set))
            out.append(
                SearchResult(
                    title=anchor.name,
                    url=anchor.doi or "offline-corpus",
                    snippet=f"{anchor.citation} ({bodies})",
                )
            )
    return out


def run_one(key: str) -> dict:
    rep = REPRESENTATIVES[key]
    label, anchor_mu, x0, jacobi, period, hc = (
        rep["label"],
        rep["anchor_mu"],
        rep["x0"],
        rep["jacobi"],
        rep["period"],
        rep["hc"],
    )
    titan = cr3bp_system("Saturn", "Titan")
    c_l1_titan = _c_l1(titan.mu)

    _print_and_append(
        [
            f"=== {label} :: anchor mu={anchor_mu:.6g} ({rep['rrt_anchor_note']}) "
            f"-> Titan mu={titan.mu:.6g} (C_L1(Titan)={c_l1_titan:.6f}) ==="
        ]
    )

    t0 = time.time()
    landed = mu_step_to_system_tracking_c_l1(
        anchor_mu,
        titan,
        x0,
        jacobi,
        period,
        hc=hc,
        sign=-1.0,
        n_steps=200,
        c_margin=rep["c_margin"],
        tol=1e-10,
    )
    elapsed = time.time() - t0

    if landed is None:
        _print_and_append(
            [f"{label}: continuation to Titan mu FAILED after {elapsed:.1f}s -- clean negative."]
        )
        return {"label": label, "continued": False, "elapsed_s": elapsed}

    state0 = np.array([landed.x0, 0.0, 0.0, 0.0, landed.ydot0, 0.0])
    topo = winding_topology(titan.mu, state0, landed.period)
    nu, _lam = cp.barden_stability(titan, landed, rtol=1e-13, atol=1e-13)

    po = cp.PeriodicOrbit(
        state0=state0,
        period=landed.period,
        jacobi=landed.jacobi,
        converged=True,
        closure_residual=landed.crossing_residual,
    )
    cc_ok, cc_dj = cp.crosscheck_periodic(titan, po, closure_tol=1e-6, jacobi_tol=1e-8)

    _print_and_append(
        [
            f"{label}: landed at Titan mu in {elapsed:.1f}s. "
            f"x0={landed.x0:.9f} C={landed.jacobi:.7f} T={landed.period:.6f} TU "
            f"({landed.period * titan.t_s / 86400.0:.3f} d)",
            f"{label}: topology=({topo.k1},{topo.k2}) prograde={topo.prograde} "
            f"reaches_secondary={topo.reaches_secondary} (target ({rep['k1']},{rep['k2']}))",
            f"{label}: Barden nu={nu:.6f} stable={abs(nu) < 1.0} "
            f"crosscheck_ok={cc_ok} crosscheck_dj={cc_dj:.2e}",
        ]
    )

    result: dict = {
        "label": label,
        "continued": True,
        "elapsed_s": elapsed,
        "x0": landed.x0,
        "jacobi": landed.jacobi,
        "period": landed.period,
        "topology": (topo.k1, topo.k2),
        "prograde": topo.prograde,
        "reaches_secondary": topo.reaches_secondary,
        "nu_direct": nu,
        "stable_direct": abs(nu) < 1.0,
        "crosscheck_ok": cc_ok,
    }

    # If the landed member isn't stable/on-topology at its landed C, trace the
    # C-family at Titan mu the way #494's own Phase-3 did for Pluto-Charon
    # (sweep_31 / scan_c_family_at_mu): the family's stable window need not sit
    # exactly at the value a fixed-Jacobi walk happens to land on.
    if not (topo.k1 == rep["k1"] and topo.k2 == rep["k2"] and abs(nu) < 1.0):
        _print_and_append([f"{label}: landed member off-target/unstable -- C-sweep at Titan mu."])
        members = scan_c_family_at_mu(
            titan.mu,
            landed.x0,
            landed.jacobi,
            landed.period,
            half_crossings=hc,
            ydot0_sign=-1.0,
            dc=0.001,
            n_each=40,
        )
        on_target_stable = [
            m
            for m in members
            if m.stable
            and winding_topology(
                titan.mu,
                np.array([m.x0, 0.0, 0.0, 0.0, m.ydot0, 0.0]),
                m.period,
            ).k1
            == rep["k1"]
            and winding_topology(
                titan.mu,
                np.array([m.x0, 0.0, 0.0, 0.0, m.ydot0, 0.0]),
                m.period,
            ).k2
            == rep["k2"]
        ]
        _print_and_append(
            [
                f"{label}: C-sweep found {len(members)} members, "
                f"{len(on_target_stable)} on-target-topology AND stable."
            ]
        )
        if on_target_stable:
            best = min(on_target_stable, key=lambda m: abs(m.nu))
            landed = cp.SymmetricOrbit(
                x0=best.x0,
                ydot0=best.ydot0,
                jacobi=best.jacobi,
                t_half=best.period / 2.0,
                period=best.period,
                converged=True,
                crossing_residual=best.crossing_residual,
                n_iter=0,
            )
            state0 = np.array([landed.x0, 0.0, 0.0, 0.0, landed.ydot0, 0.0])
            topo = winding_topology(titan.mu, state0, landed.period)
            nu = best.nu
            result.update(
                x0=landed.x0,
                jacobi=landed.jacobi,
                period=landed.period,
                topology=(topo.k1, topo.k2),
                prograde=topo.prograde,
                reaches_secondary=topo.reaches_secondary,
                nu_direct=nu,
                stable_direct=abs(nu) < 1.0,
                crosscheck_ok=best.radau_djacobi < 1e-7,
                found_via_c_sweep=True,
            )
            _print_and_append(
                [
                    f"{label}: C-sweep-selected member x0={landed.x0:.9f} "
                    f"C={landed.jacobi:.7f} nu={nu:.6f} "
                    f"topology=({topo.k1},{topo.k2}) reaches_secondary={topo.reaches_secondary}"
                ]
            )

    # --- Gate (a): stability ---
    result["gate_stability_pass"] = bool(
        abs(result["nu_direct"]) < 1.0
        and result["topology"] == (rep["k1"], rep["k2"])
        and result["reaches_secondary"]
    )

    # --- Gate (b): encounter-relevance (perimoon-passage geometry) ---
    if result["gate_stability_pass"]:
        titan_radius_km = SATELLITES["Titan"].radius_eq_km
        passage = find_perimoon_passage(
            titan, state0, landed.period, titan_radius_km, rtol=1e-12, atol=1e-12
        )
        verdict = flyby_is_useful("Titan", passage.speed_rel_kms)
        result["perimoon"] = {
            "altitude_km": passage.altitude_km,
            "speed_rel_kms": passage.speed_rel_kms,
            "below_surface": passage.below_surface,
            "max_bend_deg": verdict.max_bend_deg,
            "is_useful_bend": verdict.is_useful,
            "verdict_notes": verdict.notes,
        }
        _print_and_append(
            [
                f"{label}: perimoon passage -- altitude={passage.altitude_km:.1f} km "
                f"(Titan R={titan_radius_km} km), rel. speed={passage.speed_rel_kms:.4f} km/s, "
                f"below_surface={passage.below_surface}",
                f"{label}: physical-sanity bend gate -- max_bend={verdict.max_bend_deg:.3f} deg, "
                f"is_useful(>=5deg floor)={verdict.is_useful}",
            ]
        )
        result["gate_encounter_pass"] = bool((not passage.below_surface) and verdict.is_useful)
    else:
        result["gate_encounter_pass"] = False
        _print_and_append([f"{label}: stability gate FAILED -- skipping perimoon-geometry gate."])

    # --- Literature check (offline backend, for the record; live-WebSearch
    # results documented separately in this task's notes) ---
    sig = CandidateSignature(
        primary="Saturn",
        sequence=("Titan",),
        topology_label=frozenset({"repeated-moon"}),
        vinf_per_encounter_kms=(
            (result["perimoon"]["speed_rel_kms"],) if "perimoon" in result else ()
        ),
    )
    lit = check_literature(sig, search=offline_corpus_search)
    result["lit_check_offline"] = {
        "status": lit.status,
        "citation": lit.citation,
        "confidence": lit.confidence,
    }
    _print_and_append(
        [
            f"{label}: offline literature-check status={lit.status!r} "
            f"citation={lit.citation!r} confidence={lit.confidence}"
        ]
    )

    return result


def main() -> None:
    preflight_search(
        task_no=627,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=2,
        override_reason=(
            "PILOT ONLY (per task #627 scope): exactly 2 fixed (k1,k2) representatives "
            "continued to ONE target system (Saturn-Titan), reusing #494's already-validated "
            "binary-cycler harness (positive controls in "
            "tests/search/test_ross_rt_2026_mu_family.py and "
            "tests/search/test_627_titan_pilot.py) plus one new, small, deterministic "
            "continuation function (mu_step_to_system_tracking_c_l1) -- not an open-ended "
            "discovery grid needing a timing pilot. A full multi-system sweep is explicitly "
            "OUT OF SCOPE for this task and is a follow-up recommendation only."
        ),
    )

    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(f"Task #627 Ross-RT mu-continuation Titan pilot -- run {stamp}\n")

    titan = cr3bp_system("Saturn", "Titan")
    _print_and_append(
        [
            f"Saturn-Titan CR3BP system (sourced from cr3bp_system registry): "
            f"mu={titan.mu:.10f} l_km={titan.l_km} t_s={titan.t_s:.4f}",
        ]
    )

    results = {key: run_one(key) for key in REPRESENTATIVES}

    _print_and_append(["", "=== SUMMARY ==="])
    for r in results.values():
        _print_and_append([f"{r['label']}: {r}"])

    print(f"\nAppended full log to {OUT_PATH}")


if __name__ == "__main__":
    main()
