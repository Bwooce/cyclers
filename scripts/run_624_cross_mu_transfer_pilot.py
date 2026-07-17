#!/usr/bin/env python3
"""#624: cross-mu transfer pilot for #608's generative ML seed model.

`#608` trained a linear-PCA + k-means-partitioned-Gaussian generative model
(``cyclerfinder.ml.orbit_generative.ClusteredGaussianLatentModel``) on this
project's own accumulated ``#210`` outcome-log corpus -- ALL of it Earth-Moon
CR3BP (mu=0.01215) -- and found it lifts the existing corrector's physically-
sane convergence rate ~12.25x over uniform seeding (49% vs 4%), independently
reproduced bit-for-bit. `#614` then tested two ways to improve the model
(corpus family-tagging, a nonlinear autoencoder) and found neither helps at
this corpus size -- a clean negative, but still entirely IN-DISTRIBUTION
(same mu the model was trained on).

The one open question behind `#542` (the learned-seed generative warm-start
idea): does #608's advantage TRANSFER to a mass ratio mu the model never saw
during training? This script answers that -- decisively, whichever way it
goes -- by evaluating the ALREADY-TRAINED #608 model (re-derived
deterministically from the SAME corpus/split/seed, NOT retrained, NOT a new
model/hyperparameter search) at two out-of-distribution mu values:

  * mu=0.001            -- has a sourced golden positive-control anchor
                            (Ross & Roberts-Tsoukkas 2026, Table I, Rep 1;
                            `data/golden/ross_rt_2026_cycler_families.yaml`)
  * Sun-Earth mu~3.0e-6  -- ``cr3bp_system("Sun", "Earth")``'s registered mu

**Does the model's output need rescaling for the new mu? Reasoned answer: NO,
apply it as-is.** CR3BP's own nondimensionalization (mass unit = total system
mass, length unit = primary separation, time unit = 1/mean motion) is
mu-INVARIANT by construction -- unlike, say, converting km at one system's
physical scale to km at another's, there is no unit-conversion factor between
two CR3BP mu values: (x0,...,vz0,period) are already dimensionless numbers in
the SAME normalized frame for any mu. This is also the established precedent
in this exact codebase: `#494` Phase 3 (`tests/search/test_ross_rt_2026_mu_
family.py::test_494_phase3_pluto_charon_32_branch_and_crosschecks`) seeds the
corrector at mu=0.10851 directly from the mu=0.1 Table-I anchor IC with NO
rescaling, just re-correction -- exactly the operation this script performs,
just across a much larger mu gap (0.01215 -> 0.001 is ~12x; -> 3.0e-6 is
~4000x). Feeding the model's raw decoded (state0, period) guess straight into
`correct_periodic` at the new mu is therefore the DIRECT, honest transfer
test: inventing an ad-hoc rescale would test a different, uninvestigated
method, not whether THIS model's learned density transfers.

One consequence of applying the SAME (mu-agnostic) `is_physically_sane`
Jacobi/period/z-bounds #608 used at Earth-Moon scale: those bounds are, in an
absolute sense, most tightly calibrated to Earth-Moon-scale families. This
script uses them UNCHANGED at the new mu values (per this task's explicit
"reuse #608's own evaluation protocol exactly" instruction) rather than
re-deriving mu-specific bounds -- documented here as a judgment call, not a
silent choice. Verified NOT to trivially break the test: the sourced Ross-RT
mu=0.001 Rep-1 anchor itself (C=3.032, T=14.77) falls inside these bounds, so
the yardstick is not vacuous at this mu. Since the SAME bounds are applied to
both the generated and the uniform-baseline arms at each mu, any residual
calibration mismatch is common-mode and does not bias the generated-vs-
baseline RATIO this task is decisive on.

No held-out reconstruction check (unlike #608): held-out examples are
Earth-Moon orbits, so round-tripping them through the encoder/decoder tests
same-domain reconstruction (already established by #608), not cross-mu
transfer -- out of scope here.

No catalogue writeback, no literature-novelty check -- capability transfer
evaluation only, per this task's explicit scope (`#624`, `#623` shortlist B1).
"""

from __future__ import annotations

import glob
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.ml.orbit_generative import (
    FEATURE_NAMES,
    assemble_corpus,
    fit_clustered_gaussian,
    is_physically_sane,
    uniform_bounding_box_sample,
)
from cyclerfinder.search.cr3bp_periodic import correct_periodic, ydot0_from_jacobi

_REPO_ROOT = Path(__file__).resolve().parent.parent
_OUTCOME_LOG_DIR = _REPO_ROOT / "out" / "outcome_log"
_OUT_DIR = _REPO_ROOT / "data" / "found" / "624_cross_mu_transfer_pilot"

N_GENERATE = 100
SEED = 608  # SAME seed as #608, so the re-derived model is bit-for-bit
# identical to #608's -- this is "reuse the already-trained model", not a
# new fit: fit_clustered_gaussian is a deterministic function of
# (corpus, split, seed), and the coordinating session already independently
# verified #608's numbers reproduce bit-for-bit under this exact seed.

# Ross & Roberts-Tsoukkas 2026 (arXiv:2606.29189v1), Table I, Rep 1
# (mu=0.001, k=(1,1)) -- golden positive-control anchor.
# Source: data/golden/ross_rt_2026_cycler_families.yaml
_ROSS_RT_MU = 0.001
_ROSS_RT_X0 = -0.647047499999966
_ROSS_RT_C = 3.031605708907296
_ROSS_RT_T = 14.774502790974823
_ROSS_RT_YDOT0_SIGN = -1.0  # tests/search/test_ross_rt_2026_mu_family.py _REP_SETTINGS[0]


def _corpus_paths() -> list[Path]:
    patterns = ["search_campaign*.jsonl", "jpl_corpus.jsonl"]
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(sorted(Path(p) for p in glob.glob(str(_OUTCOME_LOG_DIR / pattern))))
    if not paths:
        raise SystemExit(
            f"no outcome-log files found under {_OUTCOME_LOG_DIR} -- this script "
            "expects the project's own gitignored out/outcome_log/ corpus to be "
            "present locally (it is not committed)."
        )
    return paths


def _refine_batch(system: cr3bp.CR3BPSystem, rows: np.ndarray, *, tag: str) -> list[dict]:
    """Run each row's (state0, period) through the existing corrector (SAME as #608/#614)."""
    results = []
    for i, row in enumerate(rows):
        state0_guess = row[:6]
        period_guess = float(row[6])
        try:
            if period_guess <= 1e-6:
                raise ValueError("non-positive period guess")
            orbit = correct_periodic(system, state0_guess, period_guess, tol=1e-10, max_iter=30)
        except Exception as exc:  # a bad decoded guess must never crash the sweep
            results.append({"tag": tag, "index": i, "converged": False, "error": str(exc)})
            continue
        results.append(
            {
                "tag": tag,
                "index": i,
                "converged": bool(orbit.converged),
                "residual": orbit.closure_residual,
                "period": orbit.period,
                "jacobi": orbit.jacobi,
                "state0": orbit.state0.tolist(),
            }
        )
    return results


def _physically_sane_rate(results: list[dict]) -> tuple[float, int]:
    """Fraction of ALL rows that both converged AND landed on a physically-sane
    orbit -- SAME yardstick (`is_physically_sane`, default Earth-Moon-derived
    bounds) as #608/#614; see this script's module docstring for the
    documented judgment call on reusing it unchanged at a new mu.
    """
    n_sane = 0
    for r in results:
        if not r.get("converged"):
            continue
        if r.get("state0") is None or r.get("period") is None or r.get("jacobi") is None:
            continue
        if is_physically_sane(r["state0"], r["period"], r["jacobi"]):
            n_sane += 1
    return (n_sane / len(results) if results else 0.0), n_sane


def _ross_rt_corrector_infra_check() -> dict:
    """Positive control, part 1: confirm the GENERAL shooting corrector
    (`correct_periodic`, the same one used everywhere in this script) itself
    genuinely works at mu=0.001 by re-deriving the sourced Ross-RT Table-I
    Rep-1 anchor from its own published IC (a near-exact seed). This isolates
    "does the corrector/solver machinery work at this mu" from "does the
    MODEL's seed find this basin" -- if this check fails, a low generated/
    baseline convergence rate at mu=0.001 would be an infrastructure
    artifact, not evidence about the model.
    """
    system = cr3bp.CR3BPSystem(mu=_ROSS_RT_MU, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    ydot0 = ydot0_from_jacobi(_ROSS_RT_X0, _ROSS_RT_C, system.mu, sign=_ROSS_RT_YDOT0_SIGN)
    state0 = np.array([_ROSS_RT_X0, 0.0, 0.0, 0.0, ydot0, 0.0])
    orbit = correct_periodic(system, state0, _ROSS_RT_T, tol=1e-10, max_iter=30)
    return {
        "converged": bool(orbit.converged),
        "residual": float(orbit.closure_residual),
        "x0_err": float(abs(orbit.state0[0] - _ROSS_RT_X0)),
        "period_err": float(abs(orbit.period - _ROSS_RT_T)),
        "jacobi_err": float(abs(orbit.jacobi - _ROSS_RT_C)),
    }


def _closest_to_ross_rt_anchor(results: list[dict]) -> dict | None:
    """Positive control, part 2: among CONVERGED refined candidates (generated
    OR baseline, at mu=0.001), find whichever lands closest to the Ross-RT
    Rep-1 anchor in Jacobi constant (the axis that most determines which
    real family's basin a state sits in) and report its component-wise gaps.
    A small gap here is meaningful evidence a candidate landed in the
    ballpark of a KNOWN real family; a large gap for every candidate is not
    by itself damning (this one anchor is a single family among many that
    exist at mu=0.001), but combined with the infra check above it tells us
    whether the SEEDING (not the solver) is what's mu-limited.
    """
    best = None
    best_jacobi_gap = float("inf")
    for r in results:
        if not r.get("converged"):
            continue
        jacobi = r.get("jacobi")
        state0 = r.get("state0")
        period = r.get("period")
        if jacobi is None or state0 is None or period is None:
            continue
        gap = abs(jacobi - _ROSS_RT_C)
        if gap < best_jacobi_gap:
            best_jacobi_gap = gap
            best = {
                "tag": r["tag"],
                "jacobi_gap": gap,
                "period_gap": abs(period - _ROSS_RT_T),
                "x0_gap": abs(state0[0] - _ROSS_RT_X0),
            }
    return best


def _evaluate_mu(
    *,
    label: str,
    system: cr3bp.CR3BPSystem,
    model,
    train: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    """Generate-then-refine vs. uniform-baseline-then-refine, both at `system`'s
    mu, using the SAME already-trained #608 model and the SAME uniform-
    bounding-box baseline construction #608/#614 used (drawn from the
    Earth-Moon TRAIN feature bounding box -- the baseline's domain is
    deliberately unchanged across mu targets, since the comparison this task
    needs is "does the LEARNED density beat blind seeding at this mu", not
    "what if the baseline were also mu-aware").
    """
    print(f"\n=== mu target: {label} (mu={system.mu:.10e}) ===")

    generated = model.sample(N_GENERATE, rng=rng)
    gen_results = _refine_batch(system, generated, tag=f"{label}_generated")
    gen_rate = sum(r["converged"] for r in gen_results) / N_GENERATE
    gen_sane_rate, gen_n_sane = _physically_sane_rate(gen_results)
    print(
        f"  [generated] n={N_GENERATE}: raw convergence {gen_rate:.1%}; "
        f"converged+physically-sane {gen_sane_rate:.1%} ({gen_n_sane}/{N_GENERATE})"
    )

    baseline = uniform_bounding_box_sample(train, N_GENERATE, rng=rng)
    base_results = _refine_batch(system, baseline, tag=f"{label}_baseline")
    base_rate = sum(r["converged"] for r in base_results) / N_GENERATE
    base_sane_rate, base_n_sane = _physically_sane_rate(base_results)
    print(
        f"  [baseline]  n={N_GENERATE}: raw convergence {base_rate:.1%}; "
        f"converged+physically-sane {base_sane_rate:.1%} ({base_n_sane}/{N_GENERATE})"
    )

    if base_sane_rate > 0:
        ratio: float | None = gen_sane_rate / base_sane_rate
        print(f"  ratio generated/baseline (physically-sane): {ratio:.2f}x")
    elif gen_sane_rate > 0:
        ratio = math.inf
        print("  baseline: 0 physically-sane; generated: >0 physically-sane (infinite ratio)")
    else:
        ratio = None
        print("  BOTH generated and baseline: 0 physically-sane -- no signal either way")

    return {
        "label": label,
        "mu": system.mu,
        "generated": {
            "n": N_GENERATE,
            "raw_convergence_rate": gen_rate,
            "converged_and_physically_sane_rate": gen_sane_rate,
            "n_converged_and_physically_sane": gen_n_sane,
        },
        "baseline_uniform": {
            "n": N_GENERATE,
            "raw_convergence_rate": base_rate,
            "converged_and_physically_sane_rate": base_sane_rate,
            "n_converged_and_physically_sane": base_n_sane,
        },
        "ratio_generated_over_baseline": ratio,
        "_raw_results": gen_results + base_results,
    }


def main() -> None:
    t_start = time.time()
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] #624 cross-mu transfer pilot starting")

    paths = _corpus_paths()
    corpus = assemble_corpus(paths)  # Earth-Moon, SAME filters #608 used
    print(
        f"  scanned {corpus.n_scanned} raw lines -> {corpus.n_converged_prefilter} "
        f"converged -> {len(corpus)} unique, physically-sane "
        f"({corpus.primary}-{corpus.secondary}, mu={corpus.mu:.10f})"
    )

    rng = np.random.default_rng(SEED)
    n = len(corpus)
    idx = rng.permutation(n)
    n_train = int(n * 0.8)
    train_idx = idx[:n_train]
    train = corpus.features[train_idx]
    print(f"  split: {len(train)} train (SAME seed/split as #608 -> identical train set)")

    model = fit_clustered_gaussian(train, n_latent=5, n_clusters=8, seed=SEED)
    print(
        "  re-derived #608's PCA(5)+kmeans(8)-Gaussian model (SAME seed -> bit-for-bit identical)"
    )

    # --- Ross-RT mu=0.001 positive control, part 1: corrector infra check ---
    print("\n=== Ross-RT mu=0.001 Rep-1 golden anchor: corrector infra check ===")
    infra_check = _ross_rt_corrector_infra_check()
    print(
        f"  seeded from the published anchor IC: converged={infra_check['converged']}, "
        f"residual={infra_check['residual']:.2e}, "
        f"|x0 err|={infra_check['x0_err']:.2e}, |T err|={infra_check['period_err']:.2e}, "
        f"|C err|={infra_check['jacobi_err']:.2e}"
    )

    # --- Cross-mu evaluations ---
    mu001_system = cr3bp.CR3BPSystem(mu=0.001, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    sun_earth_system = cr3bp_system("Sun", "Earth")

    result_mu001 = _evaluate_mu(
        label="mu=0.001", system=mu001_system, model=model, train=train, rng=rng
    )
    result_sun_earth = _evaluate_mu(
        label="Sun-Earth", system=sun_earth_system, model=model, train=train, rng=rng
    )

    # --- Ross-RT positive control, part 2: closest converged candidate ---
    closest = _closest_to_ross_rt_anchor(result_mu001["_raw_results"])
    print("\n=== Ross-RT mu=0.001 Rep-1 golden anchor: closest converged candidate ===")
    if closest is not None:
        print(
            f"  closest converged candidate ({closest['tag']}): "
            f"|C gap|={closest['jacobi_gap']:.4f}, |T gap|={closest['period_gap']:.4f}, "
            f"|x0 gap|={closest['x0_gap']:.4f}  (anchor: C={_ROSS_RT_C:.4f}, "
            f"T={_ROSS_RT_T:.4f}, x0={_ROSS_RT_X0:.4f})"
        )
    else:
        print("  no converged candidate at mu=0.001 -- positive-control comparison skipped")

    elapsed = time.time() - t_start
    summary = {
        "task": "#624",
        "corpus": {
            "primary": corpus.primary,
            "secondary": corpus.secondary,
            "mu": corpus.mu,
            "n_scanned_raw_lines": corpus.n_scanned,
            "n_converged_prefilter": corpus.n_converged_prefilter,
            "n_unique_physically_sane": len(corpus),
            "n_train": len(train),
            "feature_names": list(FEATURE_NAMES),
        },
        "model": "PCA(5) + kmeans(8) per-cluster Gaussian, re-derived from #608's own seed/split "
        "(bit-for-bit identical to #608's trained model, NOT retrained)",
        "ross_rt_mu001_positive_control": {
            "anchor_source": "Ross & Roberts-Tsoukkas 2026 arXiv:2606.29189v1 Table I Rep 1",
            "anchor": {"mu": _ROSS_RT_MU, "x0": _ROSS_RT_X0, "C": _ROSS_RT_C, "T": _ROSS_RT_T},
            "corrector_infra_check": infra_check,
            "closest_converged_candidate": closest,
        },
        "mu_0_001": {k: v for k, v in result_mu001.items() if not k.startswith("_")},
        "sun_earth": {k: v for k, v in result_sun_earth.items() if not k.startswith("_")},
        "elapsed_s": elapsed,
        "seed": SEED,
    }
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")

    detail_path = _OUT_DIR / "refine_results.jsonl"
    with detail_path.open("w") as fh:
        for r in result_mu001["_raw_results"] + result_sun_earth["_raw_results"]:
            fh.write(json.dumps(r) + "\n")

    print(f"\n[{time.strftime('%Y-%m-%dT%H:%M:%S')}] done in {elapsed:.1f}s")
    print(f"  wrote {summary_path}")
    print(f"  wrote {detail_path}")

    print("\nHONEST HEADLINE (converged+physically-sane rate, the decisive metric):")
    for label, res in (("mu=0.001", result_mu001), ("Sun-Earth", result_sun_earth)):
        g = res["generated"]["converged_and_physically_sane_rate"]
        b = res["baseline_uniform"]["converged_and_physically_sane_rate"]
        ratio = res["ratio_generated_over_baseline"]
        ratio_str = (
            f"{ratio:.2f}x" if isinstance(ratio, float) and math.isfinite(ratio) else str(ratio)
        )
        print(f"  {label:10s}: generated {g:.1%} vs baseline {b:.1%}  (ratio: {ratio_str})")
    print(
        "\n  Earth-Moon in-distribution reference (#608/#614): generated 49% vs baseline 4% "
        "(~12.25x)"
    )


if __name__ == "__main__":
    main()
