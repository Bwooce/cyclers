"""#628: reusable seed-generation API built on #608's/#624's generative orbit model.

`#608` built a statistical generative seed model
(``cyclerfinder.ml.orbit_generative.ClusteredGaussianLatentModel``) trained on
this project's own accumulated `#210` corrector-outcome-log corpus at
Earth-Moon CR3BP (mu=0.01215), and measured a real in-distribution lift in
converged-AND-physically-sane rate over blind uniform seeding
(`#642`-corrected: ~13-27x, not the original stale 12.25x headline -- see
below). `#624` initially reported this lift TRANSFERRING to mass ratios the
model never trained on -- 30x at mu=0.001, 3.5x at Sun-Earth mu~3.0e-6 -- but
`#642` found BOTH cross-mu results were 100% degenerate L4/L5
Lagrange-point-equilibrium false positives (fixed points with zero velocity
that trivially satisfy any periodicity residual): 0 real generated orbits at
either mu, reproduced independently twice each. `#642`'s Opus adjudication
additionally found the collapse STRUCTURAL (the model has no
mu-conditioning), not a small-N sampling artifact a larger re-pilot could
fix, so `#643` purged both falsified anchors from ``LIFT_ANCHORS`` --
:func:`expected_lift_for_mu` now returns the honest, current contract (see
its own docstring). Both `#608` and `#624` were evaluation-only: one-off
``scripts/run_608_*.py`` / ``scripts/run_624_cross_mu_transfer_pilot.py``
scripts that re-derive the model from scratch and print results, not
something a real discovery task can ``import`` and call.

**`#642`/`#643` correction (read before trusting any lift number)**: the
original 12.25x/30x/3.5x headline numbers were measured before
``is_physically_sane`` rejected trivial Lagrange-point-equilibrium "orbits".
`#642` re-derived all three from #608's/#624's own saved raw converged states
through the now-fixed filter and found SEVERE contamination -- not a minor
correction (see the numbers above). This module's ``is_physically_sane``
import (via ``cyclerfinder.ml.orbit_generative``) now rejects equilibria by
default, so :func:`generate_and_refine_seeds` is CORRECT going forward, and
`#643` purged the two falsified cross-mu anchors from ``LIFT_ANCHORS`` so
:func:`expected_lift_for_mu` no longer returns a number derived from them.
See `#642`'s/`#643`'s ``data/OUTSTANDING.md`` bullets and
``data/found/642_equilibrium_contamination_audit/summary.json`` for the full
corrected numbers.

This module is the productionized version: a clean function
(:func:`generate_and_refine_seeds`) that, given a target mu (directly, or via
a named ``(primary, secondary)`` pair), returns N candidate
``(state0, period)`` seeds sampled from `#608`'s trained density and refined
by the SAME existing ``cr3bp_periodic.correct_periodic`` corrector `#608`/
`#624` used -- suitable for a real sweep/discovery script to call directly.

Three things this module does NOT do, deliberately:

* It does NOT retrain the model. :func:`get_default_model` re-derives
  `#608`'s EXACT trained artifact (same corpus-assembly filters, same 80/20
  split, same ``fit_clustered_gaussian(n_latent=5, n_clusters=8)`` call, same
  default seed) -- bit-for-bit identical to `#608`'s/`#624`'s model, not a
  fresh fit or hyperparameter search. `#614` already showed bigger/different
  models do not help at this corpus size; that result is not revisited here.
* It does NOT hide the `#624` caveat that the model proposes a seed
  converging to SOME real, physically-sane family at similar energy -- not
  necessarily the SPECIFIC family a caller is searching for. Every
  :class:`GeneratedSeed` carries the REFINED family's own invariants
  (Jacobi constant, period, a cheap general stability index) so a caller can
  check what it actually got against what it was hoping for, rather than
  trusting a black-box "generate seeds" call.
* It does NOT assume uniform transfer across mu. `#642` found `#624`'s two
  cross-mu anchors were both L4/L5-equilibrium false positives (true
  cross-mu lift ~0x), and its Opus adjudication found the collapse
  STRUCTURAL (the model has no mu-conditioning) rather than a sampling
  artifact -- so `#643` purged them. `#649`/`#651` found a CHEAP coordinate
  transform (:func:`cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu`)
  rescues real lift at two specific re-tested mu, and this module now
  applies it automatically at any off-training mu -- but :func:`expected_lift_for_mu`
  still returns a numeric estimate ONLY for mu at/effectively-at the
  training point (the single validated `#608` in-distribution anchor,
  ~13-27x) or one of those two `#649`-validated coordinate-fix points; an
  explicit ``beyond_validated_range=True``/``estimated_lift=None``
  "cross-mu transfer unvalidated" signal covers everything else -- never a
  numeric lift interpolated or extrapolated past validated data.

**Integration pattern (`#634` design read, 2026-07-18)**: this module is
consumed as a STANDALONE LIBRARY CALLABLE -- a discovery script imports
:func:`generate_and_refine_seeds` explicitly where its author judges blind
seeding would otherwise be used. There is deliberately NO
``--seed-source=generative`` CLI convention wired into ``scripts/run_*.py``
(see `#634`'s bullet in ``data/OUTSTANDING.md`` for the full survey; short
version: the ``run_*.py`` population is write-once-per-task with no shared
argument parser to hook a flag into, and flag-style interchangeability would
undermine caveat (b) above by inviting scripts to treat a generative seed as a
drop-in for a targeted analytic/anchored one). Discoverability is handled by a
cross-reference in ``search/cr3bp_seed_generator.py``'s module docstring --
the module a future seed-hunting author finds first.

**Corpus-circularity warning**: the default model is trained on the `#210`
outcome logs under ``out/outcome_log/`` (grown by
``scripts/search_campaign_daemon.py``), and ``correct_periodic`` -- which this
module calls per seed -- auto-logs every solve whenever
``CYCLERFINDER_OUTCOME_LOG`` is set. A caller running generative seeds with
that env var pointed at the training corpus would feed the model's own outputs
back into future retrains (a distribution feedback loop that silently narrows
the learned density). Unset the env var for generative-seeded runs, or log to
a separate shard excluded from :func:`default_corpus_paths`' glob.

See also :func:`calibrate_cluster_weights_for_mu` -- a cheap, optional,
per-target-mu RECALIBRATION (reweighting the existing trained density's
8 cluster mixture weights by their own observed success rate at the target
mu, not a retrain) explored per this task's step 4. See that function's
docstring for the honest result: a real pilot at the hardest-hit target
(Sun-Earth) made the physically-sane rate WORSE (8%->4%), not better,
most likely because a 15-draws-per-cluster probe budget is too noisy to
distinguish genuine per-cluster unsuitability from sampling luck -- a
documented negative result, not forced into a claimed improvement.

**`#649`/`#651` correction (2026-07-19): cross-mu is no longer purely
"unvalidated"; a coordinate fix now rescues real value.** `#649` found that
feeding `#628`'s raw Earth-Moon-shaped model output directly to a foreign-mu
corrector (the behavior described in bullet 3 above) is not the only option:
reinterpreting the raw ``(state0, period)`` guess through a mu-INVARIANT
coordinate (`#629`'s rho=(C-3)/(C_L1(mu)-3) scaled-energy quantity, plus an
exact Hill-radius length/velocity-direction scale, then solving for the
velocity MAGNITUDE that hits the target-mu Jacobi constant exactly --
:func:`cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu`)
BEFORE handing the guess to the corrector rescues real, reproduced lift at
the two mu targets `#642` found totally collapsed: ~13.5% vs 0% baseline at
mu=0.001, ~4.0% vs 1.0% baseline at Sun-Earth mu (both combined over two
independent rng draws -- see `#649`'s own ``data/OUTSTANDING.md`` RESULT
paragraph for the full numbers/construction/caveats). `#651` wires this
transform into :func:`generate_and_refine_seeds` itself: whenever a caller
requests a mu meaningfully different from ``TRAINING_MU`` (the SAME
``VALIDATED_DELTA_LOG10_MU`` gate `#643` already established), every raw
generated seed is coordinate-transformed before refinement, automatically --
a caller does not need to call the transform module directly. This is
NOT a claim that cross-mu transfer now works everywhere: `#649` tested
EXACTLY two mu points (mu=0.001, Sun-Earth mu~3.003e-6) and its own bullet
explicitly cautions broader mu coverage is untested -- :func:`expected_lift_for_mu`
therefore returns a validated numeric lift ONLY at those two exact points
(the transform is still applied at any other off-training mu, since it costs
nothing and cannot be worse than the raw approach's proven-0% collapse, but
no numeric lift is claimed there). `#649`'s own honest caveats apply
in production exactly as they did in its pilot: hits cluster on ~1 dominant
real family per mu (not a diverse set -- inspect each :class:`GeneratedSeed`'s
own ``jacobi``/``period`` to see this), no literature-novelty check is run
here at all (never was, cross-mu or not), and period is left UNCHANGED by
the transform as a leading-order Hill-scaling approximation, not an exact
rescale.
"""

from __future__ import annotations

import glob
import math
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import CR3BPSystem, cr3bp_system
from cyclerfinder.ml.cross_mu_coordinate_transform import transform_seed_to_target_mu
from cyclerfinder.ml.orbit_generative import (
    DEFAULT_JACOBI_BOUNDS,
    DEFAULT_MAX_ABS_VZ0,
    DEFAULT_MAX_ABS_Z0,
    DEFAULT_PERIOD_BOUNDS,
    ClusteredGaussianLatentModel,
    OrbitCorpus,
    assemble_corpus,
    fit_clustered_gaussian,
    is_physically_sane,
)
from cyclerfinder.search.cr3bp_periodic import correct_periodic

# `cr3bp` is re-exported (not just used internally, e.g. line ~377's
# ``cr3bp.propagate``) because the test suite monkeypatches
# ``seed_generation.cr3bp.propagate`` directly to exercise the propagation-
# failure path without a real corpus. `__all__` makes that re-export explicit
# for mypy's strict-mode implicit-reexport check (same pattern as
# `scripts/certify_610_proteus_bend_interval.py`).
__all__ = ["cr3bp"]

_REPO_ROOT = Path(__file__).resolve().parents[3]
_DEFAULT_OUTCOME_LOG_DIR = _REPO_ROOT / "out" / "outcome_log"

# `#608`'s Earth-Moon training mass ratio -- the model's in-distribution anchor
# and the reference point every ``delta_log10_mu`` is measured from.
TRAINING_MU: float = 0.01215
DEFAULT_MODEL_SEED: int = 608
DEFAULT_N_LATENT: int = 5
DEFAULT_N_CLUSTERS: int = 8


# ---------------------------------------------------------------------------
# Model acquisition -- reproduces #608's exact training artifact, never a
# fresh/alternative fit.
# ---------------------------------------------------------------------------


def default_corpus_paths(outcome_log_dir: Path | None = None) -> list[Path]:
    """Locate the project's own `#210` outcome-log corpus files.

    Same glob `#608`'s and `#624`'s scripts used (``search_campaign*.jsonl``,
    ``jpl_corpus.jsonl`` under ``out/outcome_log/``), factored out here so
    both those scripts and this module's :func:`get_default_model` share one
    definition instead of two copies drifting apart.
    """
    directory = outcome_log_dir if outcome_log_dir is not None else _DEFAULT_OUTCOME_LOG_DIR
    patterns = ["search_campaign*.jsonl", "jpl_corpus.jsonl"]
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(sorted(Path(p) for p in glob.glob(str(directory / pattern))))
    if not paths:
        raise FileNotFoundError(
            f"no outcome-log files found under {directory} -- get_default_model() "
            "expects this project's own gitignored out/outcome_log/ corpus to be "
            "present locally (it is not committed to the repo). Pass corpus_paths "
            "explicitly to train/evaluate against a different corpus (tests use a "
            "small synthetic one this way)."
        )
    return paths


# In-process cache: assembling the real corpus + fitting is fast (~2s for the
# current ~540k-line/54k-row corpus) but there is no reason to repeat it on
# every call within one process/sweep script, and every past task
# (#608/#614/#624) already established the deterministic-given-inputs
# contract this relies on.
_ModelCacheValue = tuple[ClusteredGaussianLatentModel, OrbitCorpus, NDArray[np.float64]]
_MODEL_CACHE: dict[tuple[tuple[str, ...], int], _ModelCacheValue] = {}


def get_default_model(
    *,
    corpus_paths: Iterable[Path] | None = None,
    seed: int = DEFAULT_MODEL_SEED,
    cache: bool = True,
) -> tuple[ClusteredGaussianLatentModel, OrbitCorpus, NDArray[np.float64]]:
    """Re-derive `#608`'s trained model, EXACTLY (not a fresh fit).

    Deterministic given ``(corpus_paths, seed)``: assembles the corpus via
    :func:`assemble_corpus`, applies the identical 80/20 rng-permutation
    train/held-out split `#608` used, and fits
    ``fit_clustered_gaussian(train, n_latent=5, n_clusters=8, seed=seed)``.
    With the default ``seed=608`` and the project's own real corpus this
    reproduces `#608`'s/`#624`'s model bit-for-bit (independently verified by
    both those tasks' own coordinating-session re-runs).

    Returns ``(model, corpus, train_features)`` -- the corpus and train split
    are returned alongside the model because several callers (novelty checks,
    the uniform-baseline comparison, recalibration probes) need the same
    train feature array the model itself was fit on.
    """
    resolved_paths = (
        tuple(corpus_paths) if corpus_paths is not None else tuple(default_corpus_paths())
    )
    key = (tuple(str(p) for p in resolved_paths), seed)
    if cache and key in _MODEL_CACHE:
        return _MODEL_CACHE[key]

    corpus = assemble_corpus(resolved_paths)
    rng = np.random.default_rng(seed)
    n = len(corpus)
    idx = rng.permutation(n)
    n_train = int(n * 0.8)
    train = corpus.features[idx[:n_train]]
    model = fit_clustered_gaussian(
        train, n_latent=DEFAULT_N_LATENT, n_clusters=DEFAULT_N_CLUSTERS, seed=seed
    )

    result = (model, corpus, train)
    if cache:
        _MODEL_CACHE[key] = result
    return result


def clear_model_cache() -> None:
    """Drop all cached models (mainly useful for tests that vary corpus_paths)."""
    _MODEL_CACHE.clear()


def resolve_system(
    *,
    mu: float | None = None,
    primary: str | None = None,
    secondary: str | None = None,
) -> CR3BPSystem:
    """Build the CR3BP system to refine seeds against.

    Either pass a named pair (``primary``/``secondary``, resolved via
    ``cyclerfinder.core.cr3bp.cr3bp_system`` -- e.g. real physical
    length/time scales for a catalogued pair), or a bare ``mu`` for a
    generic/unnamed target (mirrors `#624`'s own construction for its
    mu=0.001 target, which has no catalogued named pair): CR3BP's own
    nondimensionalization is mu-invariant, so a generic system with
    ``l_km=1.0, t_s=1.0`` is exactly as valid a refinement target as a named
    one -- only the mu value matters for the corrector.
    """
    if primary is not None or secondary is not None:
        if primary is None or secondary is None:
            raise ValueError("resolve_system: pass both primary and secondary, or neither")
        return cr3bp_system(primary, secondary)
    if mu is None:
        raise ValueError("resolve_system: must supply either mu, or (primary, secondary)")
    return CR3BPSystem(mu=mu, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)


def _system_label(system: CR3BPSystem) -> str:
    if system.primary not in ("P1", "P2") and system.secondary not in ("P1", "P2"):
        return f"{system.primary}-{system.secondary}"
    return f"mu={system.mu:.6e}"


# ---------------------------------------------------------------------------
# Lift-vs-|delta log10 mu| documentation (task step 4).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LiftAnchor:
    """One empirically measured generate-then-refine lift data point."""

    label: str
    mu: float
    lift: float
    source_task: str


# The ONLY empirical evidence this project has for how #608's lift behaves:
# #608's own in-distribution measurement at the training mu itself. #624
# originally added two cross-mu pilot points (30x at mu=0.001, 3.5x at
# Sun-Earth) claiming the lift TRANSFERS off-distribution, but `#642` found
# BOTH were built entirely from degenerate L4/L5 Lagrange-point-equilibrium
# false positives (fixed points that trivially satisfy any periodicity
# residual, not real orbits) -- 0 real generated orbits at either mu,
# independently reproduced twice each (an archival re-filter of the saved
# raw converged states, and a fresh live re-run). The true cross-mu lift is
# therefore ~0x, not 30x/3.5x. `#642`'s Opus adjudication additionally
# verified this collapse is STRUCTURAL, not a small-N sampling artifact a
# larger re-pilot could fix: the `#608` training corpus contains ZERO
# equilibria and zero mass near L4/L5 (so this is not a training-data leak),
# meaning the model -- which has no mu-conditioning -- proposes
# Earth-Moon-shaped seeds that a corrector at a foreign mu collapses onto the
# nearest large stable attractor instead. `#643` purged both falsified
# anchors accordingly. See `#642`'s/`#643`'s bullets in
# `data/OUTSTANDING.md` and `data/found/642_equilibrium_contamination_audit/
# summary.json` for the full numbers and reasoning.
#
# The surviving anchor's value (13.5, not #608's original stale 12.25) is
# `#642`'s fresh LIVE re-derivation of the productionized API (27% generated
# vs. 2% baseline physically-sane rate, N=100, independent rng draws from
# #608's original run) -- the single most defensible point estimate `#642`
# produced. The true in-distribution range across ALL of `#642`'s
# re-derivations (an archival re-filter of #608's own saved states, plus this
# fresh live re-run) is ~13-27x -- noisier than #608 originally reported,
# because genuine (non-equilibrium) baseline hits are rare at this N (the
# archival re-filter found 0/100 real baseline hits, making that specific
# ratio undefined; the live re-run found 2/100).
LIFT_ANCHORS: tuple[LiftAnchor, ...] = (
    LiftAnchor("earth_moon_in_distribution", TRAINING_MU, 13.5, "#608"),
)

# How close a target mu must be (in |delta log10 mu|) to TRAINING_MU to be
# treated as "in-distribution" and get the validated numeric lift above,
# rather than :func:`expected_lift_for_mu`'s explicit cross-mu-unvalidated
# signal. Deliberately tiny -- effectively "the training mu itself, plus the
# real Earth-Moon system's own mu (``cr3bp_system("Earth", "Moon").mu``,
# ~0.0121506) differing negligibly (~2e-5 in delta_log10_mu) from the literal
# ``TRAINING_MU`` constant #608 fit its model against" -- not a judgment call
# about how far transfer might plausibly extend: `#642` found ZERO validated
# evidence at ANY nonzero distance from the training mu (the CLOSEST ever
# tested point, mu=0.001 at delta~1.08, is the one that collapsed to 0 real
# orbits), so there is no data-supported basis for a wider window.
VALIDATED_DELTA_LOG10_MU: float = 1e-3


def delta_log10_mu(target_mu: float, training_mu: float = TRAINING_MU) -> float:
    """|log10(target_mu) - log10(training_mu)| -- the distance axis `#624` used."""
    return abs(math.log10(target_mu) - math.log10(training_mu))


# `#649`/`#651`: TWO SPECIFIC cross-mu points where a coordinate transform
# (:func:`cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu`,
# `#649`'s rho=(C-3)/(C_L1(mu)-3) energy match + exact Hill-radius
# position/velocity-direction scale + velocity-magnitude solve) was
# INDEPENDENTLY VALIDATED -- twice, on two separate rng draws, at each point
# -- to rescue real, substantial physically-sane lift from `#642`'s confirmed
# 0% collapse of the raw (untransformed) approach. This is deliberately NOT a
# general "coordinate-fixed cross-mu lift" formula: `#649`'s own bullet in
# `data/OUTSTANDING.md` explicitly cautions that broader mu coverage beyond
# these two exact points is untested. `mu_0.001`'s ``lift`` is ``math.inf``
# (not merely large) because the fresh in-run baseline scored literally 0/200
# real hits at that mu both times -- an undefined/infinite ratio is the
# honest number, not a fabricated finite one.
CROSS_MU_COORDINATE_FIX_ANCHORS: tuple[LiftAnchor, ...] = (
    # 27/200 (13.5%) transformed-generated vs 0/200 (0.0%) fresh baseline,
    # combined over two independent rng draws (run 1: 15/100 vs 0/100; run 2,
    # --draw-seed=12649: 12/100 vs 0/100).
    LiftAnchor("mu_0.001_coordinate_fixed", 0.001, math.inf, "#649"),
    # 8/200 (4.0%) transformed-generated vs 2/200 (1.0%) fresh baseline,
    # combined over two independent rng draws (run 1: 5/100 vs 1/100; run 2:
    # 3/100 vs 1/100).
    LiftAnchor("sun_earth_coordinate_fixed", cr3bp_system("Sun", "Earth").mu, 4.0, "#649"),
)

# How close target_mu must be (in |delta log10 mu|, evaluated against that
# SPECIFIC anchor's own mu, not TRAINING_MU) to one of the two points above
# to be treated as "the #649-validated coordinate-fix point" and get a
# numeric lift estimate, rather than the honest cross-mu-unvalidated signal.
# Deliberately tiny for the same reason `VALIDATED_DELTA_LOG10_MU` is tiny:
# `#649` tested exactly these two mu values and nothing else -- this is a
# point lookup (allowing for the Sun-Earth mu's own floating-point
# derivation via `cr3bp_system` varying in the last few digits across
# equivalent call sites), not a window to interpolate or extrapolate within.
CROSS_MU_COORDINATE_FIX_TOLERANCE: float = 1e-6


@dataclass(frozen=True)
class LiftEstimate:
    target_mu: float
    delta_log10_mu: float
    # `#643`: no longer ALWAYS a number when ``beyond_validated_range`` is
    # ``True``. `#649`/`#651` narrow this further: ``estimated_lift`` is
    # non-``None`` in exactly two cases -- (1) at/near ``training_mu``
    # (``beyond_validated_range=False``), or (2) at one of the two SPECIFIC
    # `#649`-validated coordinate-fix mu points (``beyond_validated_range``
    # stays ``True`` there -- it is still not in-distribution parity -- but
    # ``coordinate_fix_anchor_label`` is set to say which one). Anywhere
    # else, both ``estimated_lift`` and ``coordinate_fix_anchor_label`` are
    # ``None`` -- see :func:`expected_lift_for_mu`.
    estimated_lift: float | None
    method: str
    beyond_validated_range: bool
    caveat: str
    # `#649`/`#651`: set to the matching ``CROSS_MU_COORDINATE_FIX_ANCHORS``
    # label (e.g. ``"mu_0.001_coordinate_fixed"``) when ``target_mu`` is one
    # of the two `#649`-validated points, else ``None``. Lets a caller
    # distinguish "validated coordinate-fix point" from "in-distribution" or
    # "genuinely unvalidated" without string-matching ``method``/``caveat``.
    coordinate_fix_anchor_label: str | None = None


def expected_lift_for_mu(target_mu: float, *, training_mu: float = TRAINING_MU) -> LiftEstimate:
    """Honest, DOCUMENTED expectation for the generate-then-refine physically-sane
    convergence lift at ``target_mu``.

    **`#643` correction**: this function previously interpolated/extrapolated
    a numeric lift through 3 empirical anchors (`#608` in-distribution
    12.25x, `#624`'s 30x at mu=0.001 and 3.5x at Sun-Earth). `#642` found
    `#624`'s two cross-mu anchors were both built from degenerate L4/L5
    Lagrange-point-equilibrium false positives -- 0 real generated orbits at
    either mu, reproduced independently twice each -- and its Opus
    adjudication found the collapse STRUCTURAL (the model has no
    mu-conditioning), not a sampling artifact a larger re-pilot could fix.
    Those two anchors are REMOVED from ``LIFT_ANCHORS``; this function no
    longer fabricates a numeric lift for mu it has no validated evidence
    about. Current behavior:

    * For ``target_mu`` at or effectively at ``training_mu``
      (``delta_log10_mu <= VALIDATED_DELTA_LOG10_MU``): returns the single
      validated `#608` in-distribution anchor (~13.5x, `#642`-corrected
      live-re-run point estimate; the true range across all of `#642`'s
      re-derivations is ~13-27x) and ``beyond_validated_range=False``. This
      branch is UNCHANGED by `#649`/`#651` (see below) -- verified by an
      explicit regression test.
    * For anything meaningfully different from the training mu (including
      `#624`'s own former mu=0.001 and Sun-Earth anchor points): returns
      ``estimated_lift=None`` and ``beyond_validated_range=True``, with a
      ``caveat`` explicitly stating cross-mu transfer is UNVALIDATED --
      in fact, the only two points ever tested there both collapsed to
      ~0x -- and that a caller should prefer the uniform baseline's own
      measured rate at the target mu over any generative-model lift claim.

    **`#649`/`#651` addition**: the statement above ("only two points ever
    tested there both collapsed to ~0x") described the RAW, untransformed
    approach `#642` tested. `#649` re-tested those SAME two exact points
    (mu=0.001, Sun-Earth mu~3.003e-6) with a coordinate transform applied
    first (`#649`'s rho-matched Hill-radius scaling --
    :func:`cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu`,
    now wired into :func:`generate_and_refine_seeds` itself for any
    off-training mu) and found REAL, reproduced lift -- ~13.5x-equivalent at
    mu=0.001 (baseline was 0%, so the ratio is technically infinite, not
    merely large) and ~4.0x at Sun-Earth, both combined over two independent
    rng draws. So there is now a THIRD case, inserted between the two above:
    for ``target_mu`` within ``CROSS_MU_COORDINATE_FIX_TOLERANCE`` of one of
    ``CROSS_MU_COORDINATE_FIX_ANCHORS``' two specific points, this function
    returns that anchor's validated numeric ``estimated_lift`` (with
    ``coordinate_fix_anchor_label`` set and ``beyond_validated_range`` still
    ``True`` -- this is NOT full in-distribution parity, just a real,
    validated rescue at that one exact point). Deliberately still NOT a
    general cross-mu formula: `#649`'s own bullet in ``data/OUTSTANDING.md``
    explicitly cautions broader mu coverage beyond these two points is
    untested, so anything else still falls through to the honest
    unvalidated final branch below, unchanged.

    See `#642`/`#643`/`#649`/`#651` in ``data/OUTSTANDING.md`` for the full
    numbers and reasoning behind this narrowing.
    """
    anchor = LIFT_ANCHORS[0]
    delta = delta_log10_mu(target_mu, training_mu)

    if delta <= VALIDATED_DELTA_LOG10_MU:
        return LiftEstimate(
            target_mu=target_mu,
            delta_log10_mu=delta,
            estimated_lift=anchor.lift,
            method=(
                f"validated in-distribution anchor ({anchor.label}, "
                f"{anchor.source_task}, #642-corrected)"
            ),
            beyond_validated_range=False,
            caveat=(
                f"At/near the training mu ({training_mu!r}): a real, substantial "
                "in-distribution lift is validated (~13-27x across #642's "
                f"re-derivations, {anchor.lift:.1f}x point estimate here). This "
                "does NOT extend to other mu -- #624's own cross-mu claim at "
                "this same lift level was falsified by #642; see the "
                "unvalidated branch this function returns for any other mu."
            ),
        )

    for cf_anchor in CROSS_MU_COORDINATE_FIX_ANCHORS:
        cf_delta = delta_log10_mu(target_mu, cf_anchor.mu)
        if cf_delta <= CROSS_MU_COORDINATE_FIX_TOLERANCE:
            lift_str = (
                "infinite (baseline scored 0%)"
                if math.isinf(cf_anchor.lift)
                else (f"{cf_anchor.lift:.1f}x")
            )
            return LiftEstimate(
                target_mu=target_mu,
                delta_log10_mu=delta,
                estimated_lift=cf_anchor.lift,
                method=(
                    f"#649-validated coordinate-transform-rescued cross-mu anchor "
                    f"({cf_anchor.label}, {cf_anchor.source_task}) -- only valid "
                    "because generate_and_refine_seeds applies "
                    "cross_mu_coordinate_transform.transform_seed_to_target_mu to "
                    "every raw seed automatically at this mu (see its own docstring)"
                ),
                beyond_validated_range=True,
                caveat=(
                    f"This exact mu is one of the TWO SPECIFIC cross-mu points #649 "
                    "independently re-tested (twice each, different rng draws) with "
                    "a Hill-radius/scaled-energy coordinate transform applied before "
                    f"refinement, rescuing real lift from #642's confirmed 0% "
                    f"collapse of the raw approach: {lift_str} transformed-generated "
                    "vs. baseline. This does NOT generalize to nearby or other mu -- "
                    "#649's own honest caveats: (1) hits cluster on ~1 dominant "
                    "family per mu, not a diverse set; (2) no literature-novelty "
                    "check is run; (3) period is left unchanged by the transform as "
                    "a leading-order Hill-scaling approximation, not exact; (4) only "
                    "these two exact mu points were re-tested -- broader coverage "
                    "is UNVALIDATED, see the fallback branch for any other mu. See "
                    "#649/#651 in data/OUTSTANDING.md for the full numbers."
                ),
                coordinate_fix_anchor_label=cf_anchor.label,
            )

    return LiftEstimate(
        target_mu=target_mu,
        delta_log10_mu=delta,
        estimated_lift=None,
        method="cross-mu transfer UNVALIDATED -- no numeric estimate returned",
        beyond_validated_range=True,
        caveat=(
            "Cross-mu transfer unvalidated -- prefer the uniform baseline's own "
            "measured rate at this mu, see #642. #624 originally claimed 30x at "
            "mu=0.001 and 3.5x at Sun-Earth, but #642 found BOTH were 100% "
            "degenerate L4/L5 Lagrange-point-equilibrium false positives (0 real "
            "generated orbits at either mu, in two independent trials each) -- "
            "the true cross-mu lift is ~0x, and #642's Opus adjudication found "
            "this collapse STRUCTURAL (the model has no mu-conditioning), not an "
            "N-too-small sampling artifact. Do not treat this model as validated "
            "for seed generation this far from its training mu without a "
            "genuinely new mu-conditioned model. #649 found a coordinate transform "
            "(see cyclerfinder.ml.cross_mu_coordinate_transform, now applied "
            "automatically by generate_and_refine_seeds at this mu) rescues real "
            "lift, but ONLY validated it at two specific mu (mu=0.001, Sun-Earth) "
            "-- this mu is not one of those, so no numeric lift is claimed here. "
            "See #642/#643/#649/#651 in data/OUTSTANDING.md for the full reasoning."
        ),
    )


# ---------------------------------------------------------------------------
# Cheap general stability proxy (task step 5: "stability if cheap to compute").
# ---------------------------------------------------------------------------


def stability_index(
    system: CR3BPSystem, state0: NDArray[np.float64], period: float
) -> tuple[float | None, str]:
    """Cheap, GENERAL stability proxy: spectral radius (max |eigenvalue|) of the
    full-period monodromy matrix.

    This is deliberately NOT ``cr3bp_periodic.barden_stability``'s specific
    ``nu`` parameter -- that function requires a perpendicular-x-axis-crossing
    SYMMETRIC orbit IC (``SymmetricOrbit``, form ``(x0,0,0,0,ydot0,0)``), a
    structure this model's general 6-component decoded state does not
    guarantee. The spectral radius of the monodromy (one extra propagation
    with the STM over the already-known period -- cheap, no new corrector
    iterations) is the general-purpose CR3BP stability indicator that applies
    to ANY periodic orbit: >1 means at least one growing/decaying
    (Lyapunov-unstable) mode; ~1 across all eigenvalues means neutrally/
    marginally stable, the typical case near a family's own stability
    boundary. Returns ``(None, reason)`` if the propagation fails for any
    reason -- a stability-index failure must never crash seed generation.
    """
    try:
        arc = cr3bp.propagate(system, np.asarray(state0, dtype=np.float64), period, with_stm=True)
        if arc.stm is None:
            return None, "monodromy unavailable (no STM returned)"
        eigs = np.linalg.eigvals(arc.stm)
        note = "max|eigenvalue| of full-period monodromy (general proxy, not Barden nu)"
        return float(np.max(np.abs(eigs))), note
    except Exception as exc:  # a stability check must never crash seed generation
        return None, f"stability computation failed: {exc}"


# ---------------------------------------------------------------------------
# Generate-then-refine: the reusable production API.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GeneratedSeed:
    """One generate-then-refine candidate.

    Carries the REFINED family's own invariants (jacobi/period/stability)
    alongside the raw seed so a caller can check what it actually got against
    what it was hoping for -- per `#624`'s honest caveat: the model proposes
    a seed that converges to SOME real, physically-sane family at similar
    energy, not necessarily the SPECIFIC family being searched for. `#649`
    found this caveat applies EVEN MORE sharply cross-mu: hits there tend to
    cluster on ~1 dominant real family per mu, not a diverse set -- inspect
    ``jacobi``/``period`` across a batch of seeds from one
    :func:`generate_and_refine_seeds` call to see this clustering directly
    (it is never hidden behind this report's aggregate rate alone).

    `#651`: ``state0_guess``/``period_guess`` are whatever was actually
    handed to ``correct_periodic`` -- the raw model draw when
    ``coordinate_transform_applied`` is ``False``, or `#649`'s
    coordinate-transformed guess (see
    :func:`cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu`)
    when it is ``True``. When the transform itself could not construct a
    valid guess for a draw (an honest ``None``, not a crash), this seed
    carries the RAW (untransformed) draw in ``state0_guess``/``period_guess``
    instead, with ``converged=False`` and ``stability_note`` explaining the
    construction failure -- counted like any other non-converged attempt,
    never silently dropped.
    """

    state0_guess: NDArray[np.float64]
    period_guess: float
    converged: bool
    physically_sane: bool
    state0: NDArray[np.float64] | None
    period: float | None
    jacobi: float | None
    residual: float | None
    stability_index: float | None
    stability_note: str
    # `#649`/`#651`: True iff this seed's target mu was far enough from
    # TRAINING_MU (per VALIDATED_DELTA_LOG10_MU) that
    # generate_and_refine_seeds applied the coordinate transform before
    # refinement -- regardless of whether the transform's own construction
    # succeeded (see the class docstring's construction-failure case).
    # Defaults to False so existing external constructors (e.g.
    # scripts/run_641_sun_jupiter_seed_census.py's own test fixtures) that
    # predate this field keep working unchanged.
    coordinate_transform_applied: bool = False


@dataclass(frozen=True)
class SeedGenerationReport:
    target_mu: float
    target_label: str
    n_requested: int
    n_attempted: int
    seeds: list[GeneratedSeed]
    n_converged: int
    n_physically_sane: int
    converged_and_physically_sane_rate: float
    target_jacobi_bounds: tuple[float, float] | None
    n_matched_target_jacobi: int | None
    lift_estimate: LiftEstimate
    # `#649`/`#651`: True iff `target_mu` was far enough from TRAINING_MU
    # (per VALIDATED_DELTA_LOG10_MU) that every seed in this report went
    # through #649's coordinate transform before refinement -- a
    # report-level convenience mirroring each GeneratedSeed's own
    # `coordinate_transform_applied` flag. Defaults to False for the same
    # external-constructor-compatibility reason as that field.
    cross_mu_coordinate_transform_used: bool = False


def generate_and_refine_seeds(
    n: int,
    *,
    mu: float | None = None,
    primary: str | None = None,
    secondary: str | None = None,
    model: ClusteredGaussianLatentModel | None = None,
    corpus_paths: Iterable[Path] | None = None,
    model_seed: int = DEFAULT_MODEL_SEED,
    target_jacobi_bounds: tuple[float, float] | None = None,
    max_oversample_factor: float = 20.0,
    jacobi_bounds: tuple[float, float] = DEFAULT_JACOBI_BOUNDS,
    period_bounds: tuple[float, float] = DEFAULT_PERIOD_BOUNDS,
    max_abs_z0: float = DEFAULT_MAX_ABS_Z0,
    max_abs_vz0: float = DEFAULT_MAX_ABS_VZ0,
    rng: np.random.Generator | None = None,
    tol: float = 1e-10,
    max_iter: int = 30,
    compute_stability: bool = True,
) -> SeedGenerationReport:
    """Generate N seeds from `#608`'s trained density and refine each with the
    EXISTING ``cr3bp_periodic.correct_periodic`` corrector at the target mu --
    the reusable, importable form of `#608`'s/`#624`'s one-off scripts.

    **No literature-novelty check is run by this function, cross-mu or not --
    never was, still isn't.** Every returned seed is an unclassified,
    physically-sane periodic orbit at best, not a vetted discovery candidate;
    a caller must run this project's own ``search/literature_check.py`` (or
    equivalent) before calling anything this function returns "novel".

    **`#649`/`#651` cross-mu behavior**: when the target mu is meaningfully
    different from ``TRAINING_MU`` (``delta_log10_mu(target_mu) >
    VALIDATED_DELTA_LOG10_MU`` -- the SAME `#643` gate
    :func:`expected_lift_for_mu` uses), every raw generated seed is passed
    through
    :func:`cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu`
    BEFORE refinement, instead of handed to the corrector raw. `#649`
    independently validated (twice, at two separate mu) that this rescues
    real, substantial physically-sane lift from a confirmed 0% collapse of
    the raw approach -- but ONLY measured it AT those two exact points
    (mu=0.001, Sun-Earth mu~3.003e-6; see ``report.lift_estimate`` for the
    validated number when ``target_mu`` is one of them, or the honest
    unvalidated signal otherwise -- :func:`expected_lift_for_mu` never
    extrapolates this to other mu). The transform is applied at ANY
    off-training mu regardless (it costs one cheap analytic solve per seed
    and cannot be worse than the raw approach's proven-0% collapse), but no
    numeric lift is CLAIMED beyond those two points. `#649`'s own honest
    caveats apply to every cross-mu call of this function: (1) hits tend to
    cluster on ~1 dominant real family per mu, not a diverse set -- inspect
    the returned seeds' own ``jacobi``/``period`` fields to see this
    directly, never hidden behind the aggregate rate alone; (2) period is
    left UNCHANGED by the transform, a leading-order Hill-scaling
    approximation, not an exact rescale. At/near ``TRAINING_MU`` this
    function's behavior is completely UNCHANGED from before `#649`/`#651`
    (verified by an explicit regression test) -- the transform is never
    applied in-distribution.

    **`#651`'s own additional honest caveat (discovered during this task's
    own regression testing, not by `#649`)**: because the coordinate
    transform succeeds at constructing genuine non-equilibrium orbits (the
    entire point of `#649`'s fix), a RARE cross-mu draw (~1 in 120 observed
    during this task's own testing) can dynamically evolve into a genuine
    close approach to the secondary body partway through propagation -- a
    case the transform's own construction check cannot catch (the initial
    condition looks unremarkable; the close approach develops later in the
    integration). ``correct_periodic``'s underlying legacy variable-step STM
    propagator (``cyclerfinder.core.cr3bp._propagate_with_stm_variable``)
    has no wall-clock or step-count budget, so such a draw can make a single
    seed's refinement take a very long time (minutes or more) instead of
    failing fast. This is a PRE-EXISTING gap in that shared propagator, not
    something `#649`'s transform or this wiring introduces -- the raw
    (untransformed) cross-mu path never triggered it only because it always
    collapsed instantly onto benign L4/L5 equilibria (per `#642`). Out of
    both `#649`'s and this task's explicit scope to fix (it would need its
    own design work: a wall-clock budget, a step-count cap, or a collision
    event on ``solve_ivp`` -- a real numerical-methods/engineering decision,
    not a mechanical one). A caller running this function at scale on
    off-training mu should apply their own external timeout/resource budget
    per call until a dedicated follow-up task hardens the propagator itself.

    Parameters
    ----------
    n :
        Number of seeds to return.
    mu, primary, secondary :
        Target system: either a bare ``mu`` (generic/unnamed pair, mirroring
        `#624`'s mu=0.001 pilot target) or a named ``(primary, secondary)``
        pair resolved via ``cr3bp_system`` -- see :func:`resolve_system`.
    model, corpus_paths, model_seed :
        Pass a pre-fit ``model`` to reuse one across many calls (e.g. sweeping
        several targets); otherwise :func:`get_default_model` re-derives
        `#608`'s exact trained artifact (cached in-process).
    target_jacobi_bounds :
        Optional caller-specified Jacobi-constant window. When given, this
        function OVERSAMPLES (drawing more than ``n`` raw candidates, up to
        ``n * max_oversample_factor`` attempts) and returns only the
        converged+physically-sane candidates whose REFINED Jacobi constant
        falls in this window, ranked by closeness to the window's midpoint.
        The returned list may be shorter than ``n`` if the budget is
        exhausted first (reported honestly via ``n_matched_target_jacobi``
        vs. ``n_requested``, never silently padded). When ``None`` (default),
        behavior is unfiltered: exactly ``n`` raw samples are drawn once and
        refined, matching `#608`'s/`#624`'s own evaluation protocol exactly
        (this is the path the regression tests pin against those tasks'
        headline numbers).
    jacobi_bounds, period_bounds, max_abs_z0, max_abs_vz0 :
        Passed to ``is_physically_sane`` -- SAME defaults `#608`/`#624` used
        (Earth-Moon-scale absolute bounds, applied unchanged at any target mu
        per `#624`'s documented judgment call: both the generated and
        baseline arms use the identical yardstick, so this is common-mode
        and does not bias a generated-vs-baseline comparison even though the
        bounds are not mu-specific).
    compute_stability :
        If ``True`` (default), compute each converged seed's general
        stability index (:func:`stability_index`) -- one extra cheap STM
        propagation per converged seed, no new corrector iterations.
    """
    system = resolve_system(mu=mu, primary=primary, secondary=secondary)
    rng = rng if rng is not None else np.random.default_rng()
    if model is None:
        model, _corpus, _train = get_default_model(corpus_paths=corpus_paths, seed=model_seed)

    # `#649`/`#651`: the SAME `#643` gate `expected_lift_for_mu` uses to
    # decide "meaningfully different from the training point" -- reused, not
    # re-invented, per this task's own scope. False (and therefore every
    # downstream branch below a no-op) at/near TRAINING_MU, which is exactly
    # how the in-distribution path stays byte-for-bit unchanged.
    apply_cross_mu_transform = delta_log10_mu(system.mu) > VALIDATED_DELTA_LOG10_MU

    def _refine_one(row: NDArray[np.float64]) -> GeneratedSeed:
        raw_state0_guess = np.asarray(row[:6], dtype=np.float64)
        raw_period_guess = float(row[6])
        state0_guess = raw_state0_guess
        period_guess = raw_period_guess
        converged = False
        sane = False
        state0 = period = jacobi = residual = idx = None
        note = "not computed (solver did not converge)"
        try:
            if apply_cross_mu_transform:
                transformed = transform_seed_to_target_mu(
                    raw_state0_guess, raw_period_guess, TRAINING_MU, system.mu
                )
                if transformed is None:
                    return GeneratedSeed(
                        state0_guess=raw_state0_guess,
                        period_guess=raw_period_guess,
                        converged=False,
                        physically_sane=False,
                        state0=None,
                        period=None,
                        jacobi=None,
                        residual=None,
                        stability_index=None,
                        stability_note=(
                            "#649 coordinate transform could not construct a valid "
                            "cross-mu seed for this draw (unrealizable velocity at "
                            "the rho-matched target Jacobi constant, or a "
                            "zero-velocity guess) -- see "
                            "cyclerfinder.ml.cross_mu_coordinate_transform."
                            "transform_seed_to_target_mu's own docstring"
                        ),
                        coordinate_transform_applied=True,
                    )
                state0_guess = transformed.state0
                period_guess = transformed.period
            if period_guess <= 1e-6:
                raise ValueError("non-positive period guess")
            orbit = correct_periodic(system, state0_guess, period_guess, tol=tol, max_iter=max_iter)
            converged = bool(orbit.converged)
            if converged:
                state0 = orbit.state0
                period = orbit.period
                jacobi = orbit.jacobi
                residual = orbit.closure_residual
                sane = is_physically_sane(
                    state0,
                    period,
                    jacobi,
                    jacobi_bounds=jacobi_bounds,
                    period_bounds=period_bounds,
                    max_abs_z0=max_abs_z0,
                    max_abs_vz0=max_abs_vz0,
                )
                if compute_stability:
                    idx, note = stability_index(system, state0, period)
                else:
                    note = "compute_stability=False"
        except Exception as exc:  # a bad decoded guess must never crash generation
            note = f"refinement failed: {exc}"
        return GeneratedSeed(
            state0_guess=state0_guess,
            period_guess=period_guess,
            converged=converged,
            physically_sane=sane,
            state0=state0,
            period=period,
            jacobi=jacobi,
            residual=residual,
            stability_index=idx,
            stability_note=note,
            coordinate_transform_applied=apply_cross_mu_transform,
        )

    all_seeds: list[GeneratedSeed] = []
    matched: list[GeneratedSeed] = []

    if target_jacobi_bounds is None:
        raw = model.sample(n, rng=rng)
        for row in raw:
            all_seeds.append(_refine_one(row))
        returned = all_seeds
        n_matched: int | None = None
    else:
        max_attempts = max(int(n * max_oversample_factor), n)
        center = 0.5 * (target_jacobi_bounds[0] + target_jacobi_bounds[1])
        attempted = 0
        batch = max(n, 10)
        while attempted < max_attempts and len(matched) < n:
            this_batch = min(batch, max_attempts - attempted)
            raw = model.sample(this_batch, rng=rng)
            for row in raw:
                seed = _refine_one(row)
                all_seeds.append(seed)
                attempted += 1
                if (
                    seed.converged
                    and seed.physically_sane
                    and seed.jacobi is not None
                    and target_jacobi_bounds[0] <= seed.jacobi <= target_jacobi_bounds[1]
                ):
                    matched.append(seed)
                if attempted >= max_attempts:
                    break
        matched.sort(key=lambda s: abs((s.jacobi or center) - center))
        returned = matched[:n]
        n_matched = len(matched)

    n_converged = sum(1 for s in all_seeds if s.converged)
    n_sane = sum(1 for s in all_seeds if s.physically_sane)
    n_attempted = len(all_seeds)

    return SeedGenerationReport(
        target_mu=system.mu,
        target_label=_system_label(system),
        n_requested=n,
        n_attempted=n_attempted,
        seeds=returned,
        n_converged=n_converged,
        n_physically_sane=n_sane,
        converged_and_physically_sane_rate=(n_sane / n_attempted if n_attempted else 0.0),
        target_jacobi_bounds=target_jacobi_bounds,
        n_matched_target_jacobi=n_matched,
        lift_estimate=expected_lift_for_mu(system.mu),
        cross_mu_coordinate_transform_used=apply_cross_mu_transform,
    )


# ---------------------------------------------------------------------------
# Optional per-target-mu recalibration (task step 4, exploratory).
# ---------------------------------------------------------------------------


def _sample_from_cluster(
    model: ClusteredGaussianLatentModel, cluster_index: int, n: int, rng: np.random.Generator
) -> NDArray[np.float64]:
    """Draw ``n`` samples from a SINGLE named cluster's Gaussian (bypassing the
    mixture-weight draw) -- used only by the recalibration probe below, which
    needs to measure each cluster's OWN success rate in isolation.
    """
    latent = rng.multivariate_normal(
        model.cluster_means[cluster_index], model.cluster_covs[cluster_index], size=n
    )
    return model.inverse_transform(latent)


@dataclass(frozen=True)
class RecalibrationResult:
    calibrated_model: ClusteredGaussianLatentModel
    per_cluster_probe_sane_rate: list[float]
    original_weights: list[float]
    calibrated_weights: list[float]
    n_probe_per_cluster: int


def calibrate_cluster_weights_for_mu(
    model: ClusteredGaussianLatentModel,
    *,
    mu: float | None = None,
    primary: str | None = None,
    secondary: str | None = None,
    n_probe_per_cluster: int = 15,
    rng: np.random.Generator | None = None,
    tol: float = 1e-10,
    max_iter: int = 30,
) -> RecalibrationResult:
    """Cheap per-target-mu RECALIBRATION: reweight (never retrain) the trained
    density's 8 cluster mixture weights by each cluster's OWN observed
    generate-then-refine physically-sane success rate AT THE TARGET MU,
    instead of `#608`'s Earth-Moon-training-corpus-derived weights (which
    reflect how common each cluster was IN THE TRAINING DATA, not how well it
    transfers to a different mu). Same means/covariances/PCA encoder as the
    input model -- only the categorical mixture weights change. This is
    exactly the "rescaling or reweighting the existing trained density's
    output" this task's step 4 asked to explore, deliberately NOT a retrain
    (`#614` already showed bigger/different models don't help at this corpus
    size).

    Probes each of the model's clusters independently (``n_probe_per_cluster``
    draws EACH, bypassing the mixture weight so every cluster gets probed
    regardless of how rarely the trained weights would otherwise sample it),
    refines with the same corrector, and Laplace-smooths the observed
    per-cluster sane-rate into a new set of mixture weights (a single
    unlucky small probe batch must not permanently zero out a cluster).

    **`#642` note**: the pilot numbers below predate `is_physically_sane`'s
    `#642` degenerate-equilibrium fix and were NOT re-derived by that audit
    (out of its explicit #608/#624 scope) -- they plausibly carry the same
    contamination pattern (`#642` found Sun-Earth generated-arm hits were
    100% equilibria in the closely-related #624 evaluation). Treat this
    result as unverified pending a fresh re-run under the fixed filter, not
    as independently confirmed.

    **`#628`'s own exploratory pilot result (real Earth-Moon-trained model,
    real Sun-Earth target, n_probe_per_cluster=15, N=100 generate-then-refine
    seeds before/after, this project's own corrector, one run -- not averaged
    over repeats, so treat as indicative not definitive): recalibration made
    it WORSE, not better.** Uncalibrated generative: 8/100 physically-sane
    (8%, vs a 2/100=2% uniform baseline measured the same run -- both in the
    same ballpark as `#624`'s own 7%-vs-2%/3.5x headline, i.e. this pilot's
    own uncalibrated arm is a sane replication). Post-recalibration: 4/100
    (4%) -- HALF the uncalibrated rate. The per-cluster probe itself
    (15 draws/cluster) found sane rates of only [0.067, 0.0, 0.0, 0.133, 0.2,
    0.2, 0.0, 0.0] across the 8 clusters -- i.e. most clusters produced ZERO
    successes in 15 draws, far too few observations to distinguish "this
    cluster is genuinely worse at Sun-Earth" from "15 draws is not enough to
    measure an 8-20% success rate reliably" (a single successful/unsuccessful
    draw swings a cluster's estimate by 1/15 ~ 6.7 percentage points). The
    Laplace-smoothed reweighting nonetheless shifted mass away from clusters
    that scored 0/15 by chance as much as by genuine unsuitability, and
    toward clusters that scored a lucky 2-3/15 -- overfitting sampling noise
    at this probe budget, not a stable structural improvement. A larger probe
    budget might resolve this (untested; `#624`'s own N=100 headline
    measurement needed 100 draws to pin down a ~3.5-8% rate confidently, so
    15/cluster is plausibly just undersized), but that is unverified, not
    assumed. **Honest verdict: recalibration as implemented here does NOT
    help at the hardest-hit target, and this task's own step-4 allowance
    ("recalibration doesn't help either" is a fine outcome) applies --
    NOT forced into a claimed improvement that the actual pilot run
    contradicts.** Left in this module as a documented, working, honestly-
    negative-result tool rather than deleted, in case a future task wants to
    revisit it with a larger probe budget.
    """
    system = resolve_system(mu=mu, primary=primary, secondary=secondary)
    rng = rng if rng is not None else np.random.default_rng()
    k = len(model.cluster_weights)
    success_counts = np.zeros(k)
    per_cluster_rate: list[float] = []

    for c in range(k):
        rows = _sample_from_cluster(model, c, n_probe_per_cluster, rng)
        n_sane = 0
        for row in rows:
            state0_guess = row[:6]
            period_guess = float(row[6])
            try:
                if period_guess <= 1e-6:
                    raise ValueError("non-positive period guess")
                orbit = correct_periodic(
                    system, state0_guess, period_guess, tol=tol, max_iter=max_iter
                )
                if orbit.converged and is_physically_sane(orbit.state0, orbit.period, orbit.jacobi):
                    n_sane += 1
            except Exception:
                pass
        success_counts[c] = n_sane
        per_cluster_rate.append(n_sane / n_probe_per_cluster)

    smoothed = (success_counts + 1.0) / (n_probe_per_cluster + 2.0)
    new_weights = smoothed / smoothed.sum() if smoothed.sum() > 0 else model.cluster_weights

    calibrated = ClusteredGaussianLatentModel(
        mean=model.mean,
        scale=model.scale,
        components=model.components,
        explained_variance_ratio=model.explained_variance_ratio,
        cluster_means=model.cluster_means,
        cluster_covs=model.cluster_covs,
        cluster_weights=np.asarray(new_weights, dtype=np.float64),
        n_latent=model.n_latent,
        seed=model.seed,
    )
    return RecalibrationResult(
        calibrated_model=calibrated,
        per_cluster_probe_sane_rate=per_cluster_rate,
        original_weights=model.cluster_weights.tolist(),
        calibrated_weights=np.asarray(new_weights, dtype=np.float64).tolist(),
        n_probe_per_cluster=n_probe_per_cluster,
    )
