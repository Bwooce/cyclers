"""Schema-v4 semantic/cross-field validation gate (spec §16.7, Task 1.3).

Covers rules that JSON Schema (Task 1.4) cannot express: cross-field
semantics and census ratchets.  Returns a list of error strings (never
raises), so callers can collect all violations at once.

Rules enforced
--------------
1. ``cycler_class in {multi-arc, non-keplerian}`` ⇒ ``orbit_elements.a_au``
   and ``orbit_elements.e`` must be absent or null (no single-ellipse
   conic element at the top-level orbit block).
2. ``cycler_class == non-keplerian`` ⇒ ``primary != "Sun"`` (CR3BP /
   rotating-frame orbits are planet-centric, not heliocentric).
3. ``period.basis`` items, when present, must each carry both ``pair``
   and ``k`` keys.
4. (v4.2, spec §16.7.9) ``trajectory.segments[].tof_days_bounds``, when
   present, must be exactly two positive numbers with ``min <= max``.  It
   is deliberately NOT required to contain ``tof_days``.
5. (v4.2, spec §16.7.9) ``source_ephemeris``, when present, must be a
   non-empty string.
6. (v4.3, spec §16.7.10) ``superseded_by`` / ``supersedes`` link targets,
   when present, must each resolve to an existing row ``id`` and must not
   point at the row's own ``id`` (referential integrity / no self-link).
   This is the cross-row check JSON Schema cannot express, and is the
   first-class consumer that justifies promoting the supersession link out
   of prose notes.

Dispatch helper
---------------
:func:`anchors_for` returns which expected-output anchor categories apply
to a given catalogue entry dict, dispatching by ``cycler_class`` per
spec §16.7.5.  Callers should use this rather than hard-coding class
checks inline.
"""

from __future__ import annotations

from typing import Any

_MULTI_CLASS = {"multi-arc", "non-keplerian"}

# Single-letter heliocentric body codes used to tell a 2-body anchor pair
# ("E-M") from a beat token ("VEM-syn") in period.pair (mirrors the loader).
_BODY_CODES: frozenset[str] = frozenset({"V", "E", "M", "S", "J"})

# V_inf physics sanity ceiling (km/s). Real Russell-Ocampo E-M cyclers reach
# 20.3 km/s, and the heliocentric max excess at the innermost body (Venus) is
# ~84 km/s retrograde; a unit error (m/s entered in a km/s field) lands at
# ~10^3. 50 km/s sits above all real catalogue data yet far below the
# 1000x unit-error class this bound exists to catch. See plan reconciliation.
_VINF_MAX_KMS = 50.0

# Tolerances for the source-independent physics identities (plan Task 2).
_A_REL_TOL = 1e-2  # a vs (peri+apo)/2, relative
_E_ABS_TOL = 5e-3  # e vs (apo-peri)/(apo+peri), absolute
_PERIOD_REL_TOL = 0.05  # period.years vs k * synodic/beat, relative
_REACH_ABS_TOL_AU = 1e-6  # numerical slack on the encounter-reach inequality


def validate_schema_invariants(rows: list[dict[str, Any]]) -> list[str]:
    """Validate schema-v4 semantic invariants across a list of raw YAML row dicts.

    Parameters
    ----------
    rows:
        List of raw row dicts (as loaded by ``yaml.safe_load`` on
        ``data/catalogue.yaml``).

    Returns
    -------
    list[str]
        All violation messages found.  Empty list means the data is clean.
        Does NOT raise — callers decide how to handle violations.
    """
    errors: list[str] = []
    # Build the id set once for the v4.3 cross-row supersession link check (Rule 6).
    known_ids = {str(r.get("id")) for r in rows if r.get("id") is not None}
    for row in rows:
        rid = row.get("id") or "<unknown>"
        cls = str(row.get("cycler_class") or "single-ellipse")
        oe = row.get("orbit_elements") or {}
        period = row.get("period") or {}

        # Rule 1: multi-arc / non-keplerian must not carry top-level a/e
        if cls in _MULTI_CLASS:
            a_au = oe.get("a_au")
            e = oe.get("e")
            if a_au is not None:
                errors.append(
                    f"{rid}: cycler_class={cls!r} must not have orbit_elements.a_au "
                    f"(got {a_au!r}); top-level a/e is only valid for single-ellipse"
                )
            if e is not None:
                errors.append(
                    f"{rid}: cycler_class={cls!r} must not have orbit_elements.e "
                    f"(got {e!r}); top-level a/e is only valid for single-ellipse"
                )

        # Rule 2: non-keplerian implies primary != "Sun"
        if cls == "non-keplerian":
            primary = str(row.get("primary") or "Sun")
            if primary == "Sun":
                errors.append(
                    f"{rid}: cycler_class=non-keplerian requires a non-Sun primary "
                    f"(CR3BP orbit is planet-centric); got primary={primary!r}"
                )

        # Rule 3: period.basis items must have pair + k
        basis = period.get("basis")
        if basis is not None:
            for i, item in enumerate(basis):
                if not isinstance(item, dict):
                    errors.append(
                        f"{rid}: period.basis[{i}] must be a dict with 'pair' and 'k', "
                        f"got {type(item).__name__!r}"
                    )
                    continue
                if "pair" not in item:
                    errors.append(f"{rid}: period.basis[{i}] is missing required key 'pair'")
                if "k" not in item:
                    errors.append(f"{rid}: period.basis[{i}] is missing required key 'k'")

        # Rule 4 (v4.2): trajectory.segments[].tof_days_bounds shape/order.
        # Note: tof_days (when present) is deliberately NOT required to lie
        # inside the bounds — e.g. the Aldrin outbound segment carries
        # tof_days=146 (circular-coplanar idealization) while Rogers et al.
        # 2012 Table 4 publishes a STOUR range of 161-172 d; both are sourced
        # framings of the same leg, so forcing containment would reject valid data.
        trajectory = row.get("trajectory") or {}
        segments = trajectory.get("segments")
        if segments is not None:
            for i, seg in enumerate(segments):
                if not isinstance(seg, dict):
                    continue
                bounds = seg.get("tof_days_bounds")
                if bounds is None:
                    continue
                seg_id = seg.get("id", i)
                if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds must be "
                        f"exactly 2 numbers [min, max] (got {bounds!r})"
                    )
                    continue
                lo, hi = bounds
                if (
                    not isinstance(lo, (int, float))
                    or not isinstance(hi, (int, float))
                    or isinstance(lo, bool)
                    or isinstance(hi, bool)
                ):
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds items must "
                        f"both be numbers (got {bounds!r})"
                    )
                    continue
                if lo <= 0 or hi <= 0:
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds values must "
                        f"both be > 0 (got {bounds!r})"
                    )
                if lo > hi:
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds must have "
                        f"min <= max (got min={lo!r}, max={hi!r})"
                    )

        # Rule 5 (v4.2): source_ephemeris, when present, must be a non-empty string.
        if "source_ephemeris" in row:
            se = row.get("source_ephemeris")
            if se is not None and (not isinstance(se, str) or se.strip() == ""):
                errors.append(
                    f"{rid}: source_ephemeris must be a non-empty string when present (got {se!r})"
                )

        # Rule 6 (v4.3): supersession link referential integrity. Targets must
        # resolve to an existing row id and must not be the row's own id.
        for field in ("superseded_by", "supersedes"):
            links = row.get(field)
            if links is None:
                continue
            if not isinstance(links, list):
                kind = type(links).__name__
                errors.append(f"{rid}: {field} must be a list of row ids (got {kind!r})")
                continue
            for target in links:
                if not isinstance(target, str) or target.strip() == "":
                    errors.append(
                        f"{rid}: {field} entries must be non-empty row-id strings (got {target!r})"
                    )
                    continue
                if target == rid:
                    errors.append(
                        f"{rid}: {field} must not reference the row's own id ({target!r})"
                    )
                    continue
                if target not in known_ids:
                    errors.append(
                        f"{rid}: {field} target {target!r} does not resolve to an existing row id"
                    )

    return errors


def _beat_years_for_token(bodies: list[str]) -> float | None:
    """Return the natural multi-body beat period (Julian years) for *bodies*,
    or ``None`` if no integer commensurability exists.

    Used to validate a beat-token ``period.pair`` (e.g. ``"VEM-syn"``). The
    body set is canonicalised by *ordering on heliocentric sma* before the
    beat search. ``multi_body_beat_days`` picks the *middle* element as the
    synodic reference, so the raw catalogue itinerary order (e.g.
    ``["E","M","V"]`` or ``["M","E","V"]``) — and even plain alphabetical
    sorting (``["E","M","V"]``) — can put a non-adjacent body in the middle
    and return no commensurability. Ordering by sma makes the radially
    middle body the reference (E for V/E/M), which is the physically correct
    choice and resolves all three VEM rows to the same 6.40-yr beat
    regardless of itinerary order (plan reconciliation, Forge R1 delta 4 /
    M8 R1 delta 3).
    """
    # Imported lazily so importing this module stays cheap (no numpy / search
    # pull-in for callers that only need the shape invariants).
    from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, PLANETS
    from cyclerfinder.search.resonance import beat_period_days, multi_body_beat_days

    unique = list(dict.fromkeys(bodies))  # de-dupe, preserve first-seen order
    if len(unique) < 2 or not all(b in PLANETS for b in unique):
        return None
    canonical = sorted(unique, key=lambda b: PLANETS[b].sma_au)
    try:
        tuples = multi_body_beat_days(canonical)
    except (ValueError, KeyError):
        return None
    if not tuples:
        return None
    return beat_period_days(canonical, tuples[0]) / DAYS_PER_JULIAN_YEAR


def validate_physical_invariants(rows: list[dict[str, Any]]) -> list[str]:
    """Validate source-INDEPENDENT physical-consistency identities (plan Task 2).

    These are golden in the strict sense: pure physics that must hold for
    every row regardless of its source, so they catch transcription / unit
    errors with no second source needed. Returns a list of violation
    strings (never raises); empty means clean.

    Checks (each applied only where the relevant fields exist / are non-null):

    * **orbit_elements** (single-ellipse a/e geometry):
      ``a ≈ (peri+apo)/2`` (rel ``1e-2``); ``e ≈ (apo-peri)/(apo+peri)``
      (abs ``5e-3``); ``0 ≤ e < 1``; ``peri ≤ a ≤ apo``.
    * **V∞** at each encounter: ``0 ≤ V∞ < 50 km/s`` (unit-error guard; see
      ``_VINF_MAX_KMS``).
    * **period**: ``years > 0`` and ``years ~ k * T`` (rel ``5e-2``) where
      ``T`` is the synodic period of a body-pair ``period.pair`` for a
      *single-ellipse* row, or the natural multi-body beat for a *beat
      token* ``period.pair`` (e.g. ``"VEM-syn"``). For a *multi-arc* row the
      synodic-integer ``k`` is an approximate label, not the true period
      (e.g. the sourced ``sanchez-net-2022-em-cycler2`` is ~8% under the
      ``k`` fit), so the strict identity is SKIPPED for that class.
    * **encounter reach** (single-ellipse with a/e + two heliocentric
      bodies): the orbit must reach each body's *encounterable* radial
      range — ``peri ≤ r_inner.aphelion`` AND ``apo ≥ r_outer.perihelion``
      — accounting for body eccentricity (Mars e≈0.093 ⇒ perihelion
      1.381 AU), not the mean sma (which would false-fail rows that
      encounter Mars near its perihelion).

    Parameters
    ----------
    rows:
        Raw row dicts (``yaml.safe_load`` on ``data/catalogue.yaml``).

    Returns
    -------
    list[str]
        All violation messages found (empty when clean).
    """
    from cyclerfinder.core.constants import PLANETS
    from cyclerfinder.search.resonance import synodic_period_years

    errors: list[str] = []
    for row in rows:
        rid = row.get("id") or "<unknown>"
        cls = str(row.get("cycler_class") or "single-ellipse")
        oe = row.get("orbit_elements") or {}
        a_au = oe.get("a_au")
        e = oe.get("e")
        peri = oe.get("perihelion_au")
        apo = oe.get("aphelion_au")

        # --- orbit-element identities (only where a/e present) -------------
        if a_au is not None and e is not None:
            if not (0.0 <= e < 1.0):
                errors.append(f"{rid}: eccentricity out of [0,1): e={e!r}")
            if peri is not None and apo is not None:
                a_mid = (peri + apo) / 2.0
                if a_mid > 0 and abs(a_au - a_mid) / a_mid > _A_REL_TOL:
                    errors.append(
                        f"{rid}: a_au={a_au!r} disagrees with (peri+apo)/2={a_mid:.4f} "
                        f"(rel tol {_A_REL_TOL})"
                    )
                denom = apo + peri
                if denom > 0:
                    e_geo = (apo - peri) / denom
                    if abs(e - e_geo) > _E_ABS_TOL:
                        errors.append(
                            f"{rid}: e={e!r} disagrees with (apo-peri)/(apo+peri)="
                            f"{e_geo:.4f} (abs tol {_E_ABS_TOL})"
                        )
                if not (peri <= a_au <= apo):
                    errors.append(
                        f"{rid}: requires peri <= a <= apo; got peri={peri!r}, "
                        f"a={a_au!r}, apo={apo!r}"
                    )

        # --- encounter reach (single-ellipse, a/e, two helio bodies) ------
        bodies = row.get("bodies") or []
        if (
            cls == "single-ellipse"
            and peri is not None
            and apo is not None
            and str(row.get("primary") or "Sun") == "Sun"
            and len(bodies) == 2
            and all(b in PLANETS for b in bodies)
        ):
            ranges = []
            for b in bodies:
                pd = PLANETS[b]
                ranges.append((b, pd.sma_au * (1 - pd.ecc), pd.sma_au * (1 + pd.ecc)))
            # inner = smaller sma, outer = larger sma
            ranges.sort(key=lambda t: t[1])
            inner_b, _inner_peri, inner_apo = ranges[0]
            outer_b, outer_peri, _outer_apo = ranges[-1]
            if peri > inner_apo + _REACH_ABS_TOL_AU:
                errors.append(
                    f"{rid}: orbit perihelion {peri!r} AU exceeds inner body "
                    f"{inner_b} aphelion {inner_apo:.4f} AU — cannot encounter {inner_b}"
                )
            if apo < outer_peri - _REACH_ABS_TOL_AU:
                errors.append(
                    f"{rid}: orbit aphelion {apo!r} AU is below outer body "
                    f"{outer_b} perihelion {outer_peri:.4f} AU — cannot encounter {outer_b}"
                )

        # --- V_inf sanity -------------------------------------------------
        for enc in row.get("vinf_kms_at_encounters") or []:
            if not isinstance(enc, dict):
                continue
            vk = enc.get("vinf_kms")
            if vk is None:
                continue
            body = enc.get("body", "?")
            if vk < 0:
                errors.append(f"{rid}: V_inf at {body} is negative ({vk!r})")
            elif vk >= _VINF_MAX_KMS:
                errors.append(
                    f"{rid}: V_inf at {body}={vk!r} km/s >= sanity ceiling "
                    f"{_VINF_MAX_KMS} km/s (suspected unit error)"
                )

        # --- period identity ---------------------------------------------
        period = row.get("period") or {}
        years = period.get("years")
        k = period.get("k")
        pair = period.get("pair")
        if years is not None:
            if years <= 0:
                errors.append(f"{rid}: period.years must be > 0; got {years!r}")
            elif k is not None and pair:
                parts = str(pair).split("-")
                is_body_pair = len(parts) == 2 and all(p in _BODY_CODES for p in parts)
                if is_body_pair:
                    if cls == "single-ellipse":
                        # Strict synodic-integer identity. Multi-arc rows carry
                        # the synodic-integer k only as an approximate label
                        # (sourced real period may differ), so skip them.
                        expected = k * synodic_period_years(parts[0], parts[1])
                        if expected > 0 and abs(years - expected) / expected > _PERIOD_REL_TOL:
                            errors.append(
                                f"{rid}: period.years={years!r} disagrees with "
                                f"k*T_syn({pair})={expected:.4f} (rel tol {_PERIOD_REL_TOL})"
                            )
                else:
                    # Beat token (e.g. "VEM-syn"): validate against the natural
                    # multi-body beat over the canonical body set.
                    beat_years = _beat_years_for_token(list(bodies))
                    if beat_years is not None and beat_years > 0:
                        expected = k * beat_years
                        if abs(years - expected) / expected > _PERIOD_REL_TOL:
                            errors.append(
                                f"{rid}: period.years={years!r} disagrees with "
                                f"k*beat({pair})={expected:.4f} (rel tol {_PERIOD_REL_TOL})"
                            )

    return errors


def validate_provenance_tags(rows: list[dict[str, Any]]) -> list[str]:
    """Validate Task-3 provenance tags WHERE PRESENT (forward-compatible).

    The per-field provenance back-fill (plan Task 3) attaches optional
    ``orbit_source`` / ``vinf_source`` / ``orbit_fidelity`` / ``vinf_fidelity``
    and an optional declared ``validation_tier`` to derivation-relevant rows.
    Those tags do not exist on the catalogue yet, so on the current data this
    function is a no-op — but it is the gate that makes them *machine-checked*
    the moment they land, so Task 3 cannot introduce a typo'd source key or a
    row that over-claims its tier.

    Checks (each only when the tag is present):

    * ``orbit_source`` / ``vinf_source`` must be :data:`SOURCE_REGISTRY` keys.
    * ``orbit_fidelity`` / ``vinf_fidelity`` must be valid :data:`Fidelity`
      tiers.
    * ``validation_tier``, when declared, must equal the tier
      :func:`classify_validation` actually computes from the row's sources +
      fidelities — a row cannot claim ``cross_validated`` while sharing a
      source or comparing across fidelities.
    """
    from cyclerfinder.data.provenance import (
        Tier,
        classify_validation,
        is_fidelity,
        is_registry_key,
    )

    errors: list[str] = []
    for row in rows:
        rid = row.get("id") or "<unknown>"
        orbit_source = row.get("orbit_source")
        vinf_source = row.get("vinf_source")
        orbit_fid = row.get("orbit_fidelity")
        vinf_fid = row.get("vinf_fidelity")
        declared = row.get("validation_tier")

        for field, val in (("orbit_source", orbit_source), ("vinf_source", vinf_source)):
            if val is not None and not is_registry_key(str(val)):
                errors.append(f"{rid}: {field}={val!r} is not a known provenance registry key")
        for field, val in (("orbit_fidelity", orbit_fid), ("vinf_fidelity", vinf_fid)):
            if val is not None and not is_fidelity(str(val)):
                errors.append(f"{rid}: {field}={val!r} is not a valid fidelity tier")

        if declared is not None:
            try:
                want = Tier(str(declared))
            except ValueError:
                errors.append(f"{rid}: validation_tier={declared!r} is not a recognised Tier value")
                continue
            same_fid = orbit_fid is not None and vinf_fid is not None and orbit_fid == vinf_fid
            got = classify_validation(
                str(orbit_source) if orbit_source is not None else None,
                str(vinf_source) if vinf_source is not None else None,
                same_fidelity=same_fid,
            )
            if got is not want:
                errors.append(
                    f"{rid}: declares validation_tier={want.value!r} but its sources/"
                    f"fidelities only support {got.value!r}"
                )

    return errors


# ---------------------------------------------------------------------------
# validation_level (schema v4.5, spec §16.7.12 / §14)
# ---------------------------------------------------------------------------

_VALIDATION_LEVELS: frozenset[str] = frozenset({"V0", "V1", "V2", "V3", "V4", "V5"})

# Mechanical-evidence registry: a row may declare a level ABOVE V0 only when the
# recorded test evidence justifies it. Each entry maps (id, level) -> the
# verbatim mechanical evidence pointer. The over-claim guard refuses any row
# claiming V1+ that is not in this registry — golden discipline: a level is only
# as high as the recorded evidence mechanically supports (when in doubt, V0).
#
# Today exactly one row earns above V0: the Aldrin outbound real-DE440 cycler
# clears spec §14 V1 (lamberthub izzo+gooding leg agreement < 1e-3 m/s AND the
# Kepler forward re-propagation residual), demonstrated with teeth by the slow
# Axis-A integration tests. Everything else (incl. the Aldrin INBOUND row, which
# no test builds/cross-checks on real ephemeris) stays V0.
_LEVEL_EVIDENCE: dict[tuple[str, str], str] = {
    ("aldrin-classic-em-k1-outbound", "V1"): (
        "spec §14 V1: real-DE440 Aldrin cycler — lamberthub izzo2015+gooding1990 "
        "per-leg agreement < V1_TOLERANCE_MPS AND Kepler forward re-propagation "
        "residual pass. tests/verify/test_agreement_lamberthub.py::"
        "test_report_includes_lamberthub_path + "
        "test_real_eph_paths_a_and_c_pass_b_flags_model_mismatch."
    ),
    ("aldrin-classic-em-k1-inbound", "V1"): (
        "spec §14 V1 (#125 Part 1): real-DE440 Aldrin INBOUND cycler — built "
        "like-for-like with the outbound twin (same loader, phase-signature "
        "resolution, and constructor), lamberthub izzo2015+gooding1990 per-leg "
        "agreement < V1_TOLERANCE_MPS AND Kepler forward re-propagation residual "
        "pass (both §14 V1 halves). Path (b) is unavailable here (single-ellipse "
        "crossing undefined for the short M->E first leg), so the agreed verdict "
        "rests on the two passing independent paths. tests/verify/"
        "test_agreement_lamberthub.py::"
        "test_inbound_real_eph_lamberthub_and_kepler_paths_pass."
    ),
    # #137 Part 1: three Russell free-return rows whose single heliocentric
    # ellipse forms a genuinely CLOSED, V_inf-continuous E->M->E cycler clear
    # spec §14 V1 mechanics LIKE-FOR-LIKE on the circular ephemeris (a
    # circular-coplanar reproduction of a circular-coplanar source): lamberthub
    # izzo2015+gooding1990 per-leg agreement < V1_TOLERANCE_MPS AND Kepler forward
    # re-propagation residual < KEPLER_REPROP_TOL_KM, with Mars-flyby V_inf
    # continuity intact (the genome honesty gate that refuses the six multi-arc
    # rows whose forced return breaks continuity by ~24 km/s). See
    # docs/notes/2026-06-07-russell12-freereturn-results.md.
    ("russell-ch4-5.30gGf3", "V1"): (
        "spec §14 V1 (#137 Part 1): circular like-for-like free-return reproduction "
        "— closed single-ellipse E->M->E arc, lamberthub izzo2015+gooding1990 per-leg "
        "agreement < V1_TOLERANCE_MPS AND Kepler forward re-propagation residual pass, "
        "Mars-flyby V_inf continuity 0.01 km/s. tests/search/"
        "test_free_return_v1_mechanics.py::"
        "test_free_return_rows_pass_section14_v1_mechanics[russell-ch4-5.30gGf3]."
    ),
    ("russell-ch4-9.94Gg3", "V1"): (
        "spec §14 V1 (#137 Part 1): circular like-for-like free-return reproduction "
        "— closed single-ellipse E->M->E arc, lamberthub izzo2015+gooding1990 per-leg "
        "agreement < V1_TOLERANCE_MPS AND Kepler forward re-propagation residual pass, "
        "Mars-flyby V_inf continuity 0.04 km/s. tests/search/"
        "test_free_return_v1_mechanics.py::"
        "test_free_return_rows_pass_section14_v1_mechanics[russell-ch4-9.94Gg3]."
    ),
    ("russell-ch4-5.75ggF3", "V1"): (
        "spec §14 V1 (#137 Part 1): circular like-for-like free-return reproduction "
        "— closed single-ellipse E->M->E arc, lamberthub izzo2015+gooding1990 per-leg "
        "agreement < V1_TOLERANCE_MPS AND Kepler forward re-propagation residual pass, "
        "Mars-flyby V_inf continuity 0.18 km/s. tests/search/"
        "test_free_return_v1_mechanics.py::"
        "test_free_return_rows_pass_section14_v1_mechanics[russell-ch4-5.75ggF3]."
    ),
    ("russell-ch4-9.353Gg2", "V1"): (
        "spec §14 V1 (#137 Part 3): deep-aphelion free-return row promoted to "
        "CLOSE-AND-MATCH by the dense phase scan (its narrow high-e t0 basin missed by "
        "the 256-point grid), then clearing §14 V1 mechanics like-for-like on the "
        "circular ephemeris — closed single-ellipse E->M->E arc, lamberthub "
        "izzo2015+gooding1990 per-leg agreement < V1_TOLERANCE_MPS AND Kepler forward "
        "re-propagation residual pass, Mars-flyby V_inf continuous. tests/search/"
        "test_free_return_v1_mechanics.py::"
        "test_free_return_rows_pass_section14_v1_mechanics[russell-ch4-9.353Gg2]."
    ),
    # 2026-06-07: the §14 V2 class-split amendment. The powered Aldrin outbound
    # clears the amended V2-POWERED gate (>=3 consecutive in-family cycles, each
    # achieving its encounters with the per-cycle maintenance applied AND bounded
    # intra-cycle drift vs the planned trajectory) — the gate the original single
    # V2 metric could not express for a per-cycle-retargeted cycler (it drifted
    # ~4.14e8 km / 3 laps under the cross-cycle rotating-frame-repeat metric, #134).
    ("aldrin-classic-em-k1-outbound", "V2"): (
        "spec §14 V2-powered (2026-06-07 class-split amendment): the powered "
        "Aldrin outbound cycler clears >=3 consecutive in-family cycles where "
        "(a) every planned encounter is achieved with the per-cycle maintenance "
        "applied (Mars-flyby V_inf continuity <= ENCOUNTER_VINF_TOL_KMS, "
        "strictly-positive in-family maintenance dV ~2.76-2.91 km/s/cycle) AND "
        "(b) intra-cycle drift vs the planned trajectory is bounded (per-leg "
        "Kepler forward-reprop residual <= INTRA_CYCLE_DRIFT_TOL_KM). The inbound "
        "twin stays V1 (its real-window solve lands on a ballistic dV~0 "
        "off-family neighbour, so 'maintenance applied' is not demonstrated). "
        "tests/verify/test_aldrin_v2_powered.py::"
        "test_aldrin_outbound_passes_v2_powered."
    ),
}


def validate_validation_level(rows: list[dict[str, Any]]) -> list[str]:
    """Validate the v4.5 ``validation_level`` tag WHERE PRESENT (spec §16.7.12).

    The level is the spec §14 gauntlet level (the highest gate a row has
    mechanically passed). This gate enforces what JSON Schema cannot:

    * the value must be one of ``V0``..``V5``;
    * a row may declare a level *above* ``V0`` only when it appears in
      :data:`_LEVEL_EVIDENCE` — i.e. recorded mechanical test evidence justifies
      it. This is the over-claim guard: no row silently promotes itself off the
      ``V0`` internal-consistency floor without a sourced, in-repo evidence
      pointer (golden discipline — when in doubt, ``V0``).

    Returns a list of violation strings (empty when clean); never raises.
    """
    errors: list[str] = []
    for row in rows:
        rid = str(row.get("id") or "<unknown>")
        level = row.get("validation_level")
        if level is None:
            continue
        level = str(level)
        if level not in _VALIDATION_LEVELS:
            errors.append(f"{rid}: validation_level={level!r} is not one of V0..V5")
            continue
        if level != "V0" and (rid, level) not in _LEVEL_EVIDENCE:
            errors.append(
                f"{rid}: validation_level={level!r} is above V0 but no recorded "
                f"mechanical evidence justifies it (not in _LEVEL_EVIDENCE) — "
                f"levels derive mechanically from recorded evidence; when in doubt, V0"
            )
    return errors


def validate_catalogue(rows: list[dict[str, Any]]) -> list[str]:
    """Run BOTH validation layers over *rows*, returning all violations.

    Single combined CI entry point (plan Task 7, closing pending #73's
    in-repo semantic half). Layers, in order:

    1. :func:`validate_schema_invariants` — cross-field *shape* / referential
       invariants JSON Schema cannot express (Rules 1-6).
    2. :func:`validate_physical_invariants` — source-independent *physics*
       identities (orbit a/e geometry, V∞ sanity, period commensurability,
       encounter reach).
    3. :func:`validate_provenance_tags` — registry-key / fidelity-tier /
       declared-tier checks on Task-3 provenance tags where present
       (no-op until the back-fill lands).
    4. :func:`validate_validation_level` — the v4.5 ``validation_level``
       over-claim guard (a row may declare V1+ only with recorded mechanical
       evidence; spec §16.7.12 / §14).

    The JSON-Schema *structural* layer (draft-2020 via the ``check-jsonschema``
    pre-commit hook) is the further, out-of-process layer; this function is the
    in-Python gate the test suite ratchets against.

    Returns the concatenation of all layers' messages (empty when clean).
    Never raises — callers decide how to surface violations.
    """
    return (
        validate_schema_invariants(rows)
        + validate_physical_invariants(rows)
        + validate_provenance_tags(rows)
        + validate_validation_level(rows)
    )


def anchors_for(entry: dict[str, Any]) -> dict[str, bool]:
    """Return which expected-output anchor categories apply to *entry*.

    Dispatches by ``cycler_class`` per spec §16.7.5.  The returned dict
    maps anchor-category names to ``True`` when that category applies and
    ``False`` when it must not be applied for the given class.

    Anchor categories
    -----------------
    ``"vinf"``
        V∞ multiset check — applicable to all classes (where
        ``vinf_kms_at_encounters`` is populated).
    ``"a_e"``
        Semi-major axis + eccentricity identity — applies only to
        ``single-ellipse``.  Must NOT be applied to multi-arc or
        non-keplerian entries (they have no single conic; the
        ``construct_resonant_cycler`` constructor must not be invoked for
        them).
    ``"period"``
        Period (years / k) check — applies to ``single-ellipse`` and
        ``multi-arc``; not meaningful for ``non-keplerian`` entries whose
        period is expressed in non-dimensional CR3BP time.
    ``"invariants"``
        Cycle-level identity (aphelion_ratio / turn_ratio /
        transit_times_days) — applies only to ``multi-arc``.
    ``"cr3bp"``
        CR3BP identity triple (jacobi_constant / period_nd /
        stability_index) — applies only to ``non-keplerian``.

    Parameters
    ----------
    entry:
        A raw YAML row dict (as loaded by ``yaml.safe_load`` on
        ``data/catalogue.yaml``).

    Returns
    -------
    dict[str, bool]
        Mapping of anchor-category names to applicability flags.
    """
    cls = str(entry.get("cycler_class") or "single-ellipse")
    return {
        "vinf": True,
        "a_e": cls == "single-ellipse",
        "period": cls in ("single-ellipse", "multi-arc"),
        "invariants": cls == "multi-arc",
        "cr3bp": cls == "non-keplerian",
    }


__all__ = [
    "anchors_for",
    "validate_catalogue",
    "validate_physical_invariants",
    "validate_provenance_tags",
    "validate_schema_invariants",
]
