"""Task #629 Phase A: SCALED-margin rerun of #627's Ross-RT mu-continuation to Titan.

#627 piloted continuing the Ross-Roberts-Tsoukkas 2026 (arXiv:2606.29189)
(k1,k2)=(1,1)/(3,3) ballistic-cycler families DOWN in mass ratio mu from their
Table-I anchors to the real Saturn-Titan mu=2.36695e-4, using
:func:`real_binary_kk_sweep.mu_step_to_system_tracking_c_l1` with a FIXED
ABSOLUTE ``c_margin`` (0.02 for (1,1), 0.005 for (3,3)) clamping C below
C_L1(mu) at every step. Both representatives came back clean negatives
(off-topology / non-convergent).

The #629 design read (docs/notes/2026-07-18-629-design-read-titan-kk-grid.md,
committed 9145d58) found strong quantitative evidence that #627's negative may
be an ARTIFACT of that fixed absolute margin, not a real physical result: the
natural energy gap to the L1 neck, C_L1(mu) - 3, shrinks as mu^(2/3) (Hill
scaling), so a margin sized for the mu=0.01215 Earth-Moon-scale corridor
(width ~0.19) is comically oversized for the mu=2.37e-4 Titan corridor (width
~0.016) -- for (1,1) the margin alone (0.02) exceeds the ENTIRE Titan corridor,
forcing the walk to C < 3 (rho < 0 -- outside the ballistic-capture band
rho in (0,1) where every known family member lives) well before Titan mu is
reached; for (3,3) the margin (0.005) equals essentially the anchor's ENTIRE
gap to its own C_L1, forcing an immediate, destabilizing C-walk at the very
first step.

This script reruns the SAME two representatives with the SAME
``mu_step_to_system_tracking_c_l1`` driver, but using its new
``c_margin_alpha`` parameter (#629 Phase A, added to
``real_binary_kk_sweep.py``): a margin SCALED to each step's own corridor
width, ``c_margin_step = alpha * (C_L1(mu_next) - 3)``, which tracks the
scaled energy ``rho = (C - 3) / (C_L1(mu) - 3)`` at a roughly constant value
(~1 - alpha) once the walk starts clamping -- matching how the two sourced
(1,1) Table-I anchors themselves sit at a near-invariant rho (~0.79-0.80)
across a 12x mu range in the design read's own analysis.

This is a DECISIVE, not exploratory, rerun -- report whichever of #629's three
possible outcomes (rescue / corridor-failure / ambiguous) actually happens.
No catalogue writeback (see the #629 OUTSTANDING.md bullet, step 7 of this
task's own instructions): any genuine rescue still needs the Phase C
literature/novelty check before any catalogue-grade claim.

Usage
-----
  uv run python scripts/run_629_phaseA_scaled_margin_titan.py --rep 11 --alpha 0.15
  uv run python scripts/run_629_phaseA_scaled_margin_titan.py --rep 33 --alpha 0.15
  uv run python scripts/run_629_phaseA_scaled_margin_titan.py --rep 33 --alpha 0.03
"""

from __future__ import annotations

import argparse
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
from cyclerfinder.search.real_binary_kk_sweep import mu_step_to_system_tracking_c_l1  # noqa: E402

OUT_PATH = (
    Path(__file__).parent.parent
    / "docs"
    / "notes"
    / "scratch"
    / "629_phase_a_scaled_margin_raw.txt"
)

_REGION_ID_TMPL = "rrt-mu-continuation-titan-phase-a-scaled-margin-{rep}-alpha{alpha}-2026-07-18"

# Same two representatives and same starting anchors #627 used (see that
# task's own script for the fold-avoidance rationale on (1,1)'s waypoint
# choice -- unchanged here, only the margin logic differs).
REPRESENTATIVES: dict[str, dict] = {
    "11": dict(
        label="(1,1)",
        anchor_mu=0.001,
        x0=-0.647047499999966,
        jacobi=3.031605708907296,
        period=14.774502790974823,
        hc=1,
        k1=1,
        k2=1,
        rrt_anchor_note=(
            "Table I Rep 1 (mu=0.001) -- the paper's own floor, same waypoint #627 used"
        ),
    ),
    "33": dict(
        label="(3,3)",
        anchor_mu=0.012150584270572,
        x0=-0.322477620583087,
        jacobi=3.183379082910527,
        period=19.503763587070285,
        hc=7,
        k1=3,
        k2=3,
        rrt_anchor_note="Table I Rep 3 (mu=0.01215) -- the task's literal starting point",
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
    """Deterministic offline literature backend (mirrors #627's own use)."""
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


def run_one(key: str, alpha: float) -> dict:
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
    c_l1_anchor = _c_l1(anchor_mu)
    rho_anchor = (jacobi - 3.0) / (c_l1_anchor - 3.0)
    rho_ceiling_predicted = 1.0 - alpha

    _print_and_append(
        [
            f"=== {label} :: alpha={alpha:.4g} :: anchor mu={anchor_mu:.6g} "
            f"({rep['rrt_anchor_note']}) -> Titan mu={titan.mu:.6g} "
            f"(C_L1(Titan)={c_l1_titan:.6f}) ===",
            f"{label}: anchor rho=(C-3)/(C_L1(anchor_mu)-3)={rho_anchor:.4f}; "
            f"predicted clamp ceiling rho=1-alpha={rho_ceiling_predicted:.4f}",
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
        c_margin_alpha=alpha,
        tol=1e-10,
    )
    elapsed = time.time() - t0

    if landed is None:
        _print_and_append(
            [
                f"{label} (alpha={alpha:.4g}): continuation to Titan mu FAILED after "
                f"{elapsed:.1f}s -- clean negative (convergence failure)."
            ]
        )
        return {
            "label": label,
            "alpha": alpha,
            "continued": False,
            "elapsed_s": elapsed,
            "rho_anchor": rho_anchor,
        }

    rho_landed = (landed.jacobi - 3.0) / (c_l1_titan - 3.0)

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
            f"{label} (alpha={alpha:.4g}): landed at Titan mu in {elapsed:.1f}s. "
            f"x0={landed.x0:.9f} C={landed.jacobi:.7f} rho={rho_landed:.4f} "
            f"T={landed.period:.6f} TU ({landed.period * titan.t_s / 86400.0:.3f} d)",
            f"{label} (alpha={alpha:.4g}): topology=({topo.k1},{topo.k2}) "
            f"prograde={topo.prograde} reaches_secondary={topo.reaches_secondary} "
            f"(target ({rep['k1']},{rep['k2']}))",
            f"{label} (alpha={alpha:.4g}): Barden nu={nu:.6f} stable={abs(nu) < 1.0} "
            f"crosscheck_ok={cc_ok} crosscheck_dj={cc_dj:.2e}",
        ]
    )

    result: dict = {
        "label": label,
        "alpha": alpha,
        "continued": True,
        "elapsed_s": elapsed,
        "rho_anchor": rho_anchor,
        "rho_predicted_ceiling": rho_ceiling_predicted,
        "x0": landed.x0,
        "jacobi": landed.jacobi,
        "rho_landed": rho_landed,
        "period": landed.period,
        "topology": (topo.k1, topo.k2),
        "prograde": topo.prograde,
        "reaches_secondary": topo.reaches_secondary,
        "nu_direct": nu,
        "stable_direct": abs(nu) < 1.0,
        "crosscheck_ok": cc_ok,
    }

    # Same escape hatch #627 used: if the landed member isn't on-topology/
    # stable at its landed C, trace the C-family at Titan mu (the family's
    # stable window need not sit exactly where the walk happened to land).
    if not (topo.k1 == rep["k1"] and topo.k2 == rep["k2"] and abs(nu) < 1.0):
        _print_and_append(
            [
                f"{label} (alpha={alpha:.4g}): landed member off-target/unstable -- "
                "C-sweep at Titan mu."
            ]
        )
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
                f"{label} (alpha={alpha:.4g}): C-sweep found {len(members)} members, "
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
                rho_landed=(landed.jacobi - 3.0) / (c_l1_titan - 3.0),
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
                    f"{label} (alpha={alpha:.4g}): C-sweep-selected member x0={landed.x0:.9f} "
                    f"C={landed.jacobi:.7f} rho={result['rho_landed']:.4f} nu={nu:.6f} "
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
                f"{label} (alpha={alpha:.4g}): perimoon passage -- altitude="
                f"{passage.altitude_km:.1f} km (Titan R={titan_radius_km} km), "
                f"rel. speed={passage.speed_rel_kms:.4f} km/s, "
                f"below_surface={passage.below_surface}",
                f"{label} (alpha={alpha:.4g}): physical-sanity bend gate -- "
                f"max_bend={verdict.max_bend_deg:.3f} deg, "
                f"is_useful(>=5deg floor)={verdict.is_useful}",
            ]
        )
        result["gate_encounter_pass"] = bool((not passage.below_surface) and verdict.is_useful)
    else:
        result["gate_encounter_pass"] = False
        _print_and_append(
            [
                f"{label} (alpha={alpha:.4g}): stability gate FAILED -- "
                "skipping perimoon-geometry gate."
            ]
        )

    # --- Literature check (offline backend, for the record) ---
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
            f"{label} (alpha={alpha:.4g}): offline literature-check status={lit.status!r} "
            f"citation={lit.citation!r} confidence={lit.confidence}"
        ]
    )

    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rep", choices=[*sorted(REPRESENTATIVES), "both"], required=True)
    parser.add_argument("--alpha", type=float, required=True)
    args = parser.parse_args()

    keys = list(REPRESENTATIVES) if args.rep == "both" else [args.rep]

    preflight_search(
        task_no=629,
        region_id=_REGION_ID_TMPL.format(rep=args.rep, alpha=args.alpha),
        method=MethodCapability(
            genome=(
                "Ross-Roberts-Tsoukkas 2026 (k1,k2)-cycler CR3BP genome (#494's fixed-Jacobi "
                "symmetric corrector + Barden stability + winding-topology classifier + "
                "independent Radau crosscheck), continued DOWN in mass ratio mu from the "
                "Earth-Moon-scale Table-I anchors toward Saturn-Titan mu=2.36695e-4"
            ),
            corrector=(
                "mu_step_to_system_tracking_c_l1 with c_margin_alpha (#629 Phase A): "
                "margin scaled to alpha*(C_L1(mu_next)-3) at every step, tracking near-"
                "constant rho=(C-3)/(C_L1(mu)-3) instead of #627's absolute c_margin"
            ),
            capability_tags=frozenset(
                {
                    "cr3bp",
                    "binary-cycler",
                    "k1k2-genome",
                    "real-moon",
                    "mu-continuation-downward",
                    "scaled-margin",
                }
            ),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=len(keys),
        override_reason=(
            "PHASE A rerun of #627's own PILOT (per task #629 scope): exactly the same 1-2 "
            "fixed (k1,k2) representatives #627 already piloted, continued to the SAME one "
            "target system (Saturn-Titan), reusing #627's already-validated harness "
            "(tests/search/test_627_titan_pilot.py) plus a small, targeted margin-scaling "
            "change (c_margin_alpha) to the same continuation function -- a decisive artifact "
            "check per the #629 design read, not an open-ended discovery grid."
        ),
    )

    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    header = (
        f"Task #629 Phase A scaled-margin rerun -- rep={args.rep} alpha={args.alpha} -- "
        f"run {stamp}\n"
    )
    if not OUT_PATH.exists():
        OUT_PATH.write_text(header)
    else:
        _append([header.rstrip("\n")])

    results = {key: run_one(key, args.alpha) for key in keys}

    _print_and_append(["", "=== SUMMARY ==="])
    for r in results.values():
        _print_and_append([f"{r['label']} alpha={args.alpha}: {r}"])

    print(f"\nAppended full log to {OUT_PATH}")


if __name__ == "__main__":
    main()
