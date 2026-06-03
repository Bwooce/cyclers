"""Catalogue loader + canonical signature + identity matcher (M7).

Spec references
---------------
* §16.1 — catalogue record schema v2 (the on-disk projection this module
  loads into :class:`CatalogueEntry`).
* §16.2 — canonical signature definition: lexicographically-minimal
  cyclic rotation of the body sequence, V∞ multiset binned to 0.05
  km/s, leg ``(a, e)`` multiset binned to 0.01 AU / 0.01, sha1 hash
  over the canonical-JSON of the in-hash fields. Pool partitioning by
  ``model_assumption`` per the section's final paragraph.
* §16.3 — three-state matcher: exact / probable-match-NEEDS-HUMAN /
  novel. Probable-match uses :func:`signature_distance` weighted L1
  below :data:`TAU_NEAR`.
* §12.2 — representation framework. Three representations
  (idealised-circular-coplanar, real-ephemeris, CR3BP) live in
  non-comparable pools; the matcher enforces this with
  :meth:`Catalog.filter` applied before signature comparison.
* §8 M7 milestone gate — "correctly tags the rediscovered E-M cyclers
  as known".

User-resolved decisions (M7 plan §8, resolved 2026-06-01)
---------------------------------------------------------
* **``sense`` field — KEEP DISTINCT.** Russell 2004 and McConaghy 2006
  catalogue outbound and inbound as separate entries; the catalogue
  has 2 Aldrin rows for this reason. The cycler dataclass now carries
  an explicit ``sense`` field (no heuristic derivation needed).
* **``period.years`` — OUT OF HASH.** Only ``k`` (integer) and ``pair``
  participate in the canonical hash. ``years = k * T_synodic(pair)``
  varies across sources (Hollister 2.13 vs spec 2.135 vs Wikipedia
  2.02 for the heliocentric ellipse); including it would silently
  break matches against literature entries that use a different
  synodic-period convention. Kept on the dataclass for display.

Plan: ``docs/phases/m7-catalogue-novelty-matching/plan.md`` §3.1.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Final, Literal

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.model.cycler import Cycler, SenseT, orbit_elements_au

# ---------------------------------------------------------------------------
# Paths + tolerances
# ---------------------------------------------------------------------------

CATALOGUE_PATH: Final[Path] = (
    Path(__file__).resolve().parent.parent.parent.parent / "data" / "catalogue.yaml"
)
"""Resolved path to ``data/catalogue.yaml``.

Matches ``tests/_catalogue_loader.py``'s pattern so the loader and the
M5 rediscovery test see the same on-disk file."""


VINF_BIN_KMS: Final[float] = 0.05
"""V∞ binning width per spec §16.2. Matches the literature reporting
precision (e.g. McConaghy 2006 5.65 ± 0.05). NOT test-tunable."""


A_BIN_AU: Final[float] = 0.01
"""Leg semi-major-axis binning width (AU) per spec §16.2. Literature
tabulates ``a`` to 2 decimal places (Rogers 2012 Table 1)."""


E_BIN: Final[float] = 0.01
"""Leg eccentricity binning width per spec §16.2."""


TAU_NEAR: Final[float] = 0.5
"""Probable-match weighted-L1 distance threshold per M7 plan §4.4.

Weighted distance accumulates ``Σ |Δvinf / 0.05| + Σ |Δ(a, e) / (0.01,
0.01)| + |Δperiod_yr / 0.5|``. A candidate below 0.5 is "one bin's
worth of mismatch accumulated" — human-review territory; above 0.5 is
a genuinely different cycler."""


# ---------------------------------------------------------------------------
# Frozen dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CanonicalSignature:
    """Spec §16.2 identity object — invariant to epoch / loop-start / noise.

    Two cyclers that compare equal under spec §16.2's identity rules
    produce signatures with bitwise-identical :attr:`hash`. Pool
    partitioning by :attr:`model_assumption` is enforced at match time
    (the field is NOT included in the hash).

    Attributes
    ----------
    bodies:
        Sorted unique tuple of body codes participating in the cycler
        (e.g. ``("E", "M")``). Body-set rotation-invariant.
    sequence_canonical:
        Lexicographically-minimal cyclic rotation of the dash-joined
        body sequence after closing-body deduplication
        (e.g. ``"E-M"``).
    sense:
        ``"outbound"`` / ``"inbound"`` / ``"n/a"`` per spec §16.2's
        "keep a separate ``sense`` field" guidance.
    period_k:
        Period in synodic multiples (integer; the exact identity
        field).
    period_years:
        Period in Julian years (informational only — NOT in hash).
        Different literature sources use slightly different synodic
        conventions (Hollister 2.13 vs spec 2.135); including this in
        the hash would silently invalidate matches across sources.
    vinf_multiset_binned:
        Sorted tuple of ``(body, vinf_binned_kms)`` over every
        encounter. Binning width :data:`VINF_BIN_KMS`. Sorted
        ascending so rotation-equivalent cyclers produce the same
        tuple.
    leg_elements_multiset_binned:
        Sorted tuple of ``(a_au_binned, e_binned)`` over every leg's
        heliocentric ellipse. Binning widths :data:`A_BIN_AU` /
        :data:`E_BIN`.
    model_assumption:
        Pool-filter tag — ``"circular-coplanar"`` /
        ``"analytic-ephemeris"`` / ``"cr3bp"``. Per spec §12.2 / §16.2
        signatures are NEVER compared across model_assumption. NOT
        included in the hash; consumed by :meth:`Catalog.filter` at
        match time.
    hash:
        ``"sha1:"`` prefix + 40-character hex digest of the
        canonical-JSON of the in-hash fields. Bitwise-stable across
        Python versions (sha1 is deterministic; the JSON serialiser
        uses sorted keys + no whitespace + ASCII-only).
    """

    bodies: tuple[str, ...]
    sequence_canonical: str
    sense: SenseT
    period_k: int
    period_years: float
    vinf_multiset_binned: tuple[tuple[str, float], ...]
    leg_elements_multiset_binned: tuple[tuple[float, float], ...]
    model_assumption: str
    hash: str


@dataclass(frozen=True)
class CatalogueEntry:
    """In-memory projection of one ``data/catalogue.yaml`` row.

    Carries every field the M7 matcher / writer / discovery workflow
    needs. Family-seed and citation-only entries set the optional
    structural fields (``vinf_kms_at_encounters``, ``legs``,
    ``orbit_elements``) to empty / None and have :attr:`signature` set
    to ``None`` — those entries remain in the catalogue but never
    participate in matching.
    """

    id: str
    name: str
    source: str
    trajectory_regime: str
    model_assumption: str
    bodies: tuple[str, ...]
    sequence_canonical: str
    sense: SenseT
    period_pair: str | None
    period_k: int | None
    period_years: float | None
    vinf_kms_at_encounters: tuple[tuple[str, float | None], ...]
    legs_tof_days: tuple[float | None, ...]
    legs_n_revs: tuple[int, ...]
    orbit_elements_a_au: float | None
    orbit_elements_e: float | None
    orbit_elements_i_deg: float | None
    primary: str
    first_published: dict[str, Any] | None
    priority_date: str | None
    our_status: str | None
    signature: CanonicalSignature | None
    signature_hash: str | None
    validation: dict[str, Any] = field(default_factory=dict)
    discovery: dict[str, Any] = field(default_factory=dict)
    discovery_run: dict[str, Any] | None = None
    raw: dict[str, Any] = field(default_factory=dict)
    """The full YAML row dict, retained so :func:`serialise_entry_yaml`
    can round-trip writeback updates without losing schema-v2 fields
    not modelled explicitly on the dataclass."""
    family: dict[str, Any] | None = None
    """Schema-v3 ``family{}`` linkage (spec §16.6.3); ``None`` =
    ungrouped. Not a signature input."""
    data_gaps: tuple[dict[str, Any], ...] = ()
    """Schema-v3 ``data_gaps[]`` known-unknown register (spec §16.6.4).
    Each dict carries ``path`` / ``kind`` / ``note`` / ``source_hint`` /
    ``todo_ref``. Enumerated by :func:`find_data_gaps`."""
    # ------------------------------------------------------------------
    # Schema v4 (spec §16.7, 2026-06-03) — additive, all optional,
    # all default to v3-compatible behaviour. NOT signature inputs.
    # ------------------------------------------------------------------
    cycler_class: str = "single-ellipse"
    """Structural kind of orbit: ``single-ellipse`` | ``multi-arc`` |
    ``non-keplerian``. Default ``single-ellipse`` = v3 behaviour."""
    orbit_elements_reference_frame: str = "heliocentric-inertial"
    """Frame tag for the orbit elements block: ``heliocentric-inertial``
    | ``planetcentric-inertial`` | ``rotating-synodic``."""
    orbit_elements_center: str = "Sun"
    """Body the orbital elements are referenced to (default ``Sun``)."""
    invariants: dict[str, Any] | None = None
    """Cycle-level identity descriptors for ``multi-arc`` orbits (e.g.
    ``aphelion_ratio``, ``transit_times_days``). ``None`` for all other
    classes."""
    cr3bp: dict[str, Any] | None = None
    """CR3BP identity tuple (``jacobi_constant``, ``period_nd``,
    ``stability_index``, …) for ``non-keplerian`` orbits. Read from
    ``orbit_elements.cr3bp`` in the YAML row. ``None`` for all other
    classes."""
    period_basis: tuple[dict[str, Any], ...] | None = None
    """Beat-period basis for n-body (VEM) cyclers: tuple of
    ``{"pair": ..., "k": ...}`` dicts from ``period.basis``. ``None``
    when not present."""
    # ------------------------------------------------------------------
    # Schema v4.1 (spec §16.7.7, 2026-06-03) — arc-type descriptors.
    # ------------------------------------------------------------------
    free_return_arcs: tuple[dict[str, Any], ...] | None = None
    """Russell's Earth-to-Earth free-return arc decomposition (spec §16.7.7).

    Distinct from ``trajectory.segments`` (encounter legs): Russell's arcs
    are Earth-to-Earth free-return trajectories; catalogue segments are
    encounter legs (E→M, M→E, E→E loop).

    Each dict carries:
    - ``arc_type``: ``"generic"`` | ``"half-rev"`` | ``"full-rev"``
    - ``resonance``: M:N string (e.g. ``"3:2"``) for full-rev; ``None`` otherwise
    - ``tof_years``: TOF in years (from descriptor first param for g/h); ``None``
      for full-rev arcs where TOF is determined by the resonance condition
    - ``raw_descriptor``: verbatim Russell token, e.g. ``"g(1.4612,526.02,Ll)"``

    ``None`` when no descriptor is available for the row (data gap, not error).
    Only meaningful for ``multi-arc`` cyclers; ``None`` on other classes.
    """

    @property
    def fully_defined(self) -> bool:
        """Derived progress flag: the orbit is completely specified.

        ``True`` iff the entry has its core physically-defining fields
        populated **and** carries no acknowledged known-unknown
        (:attr:`data_gaps`). This is a deliberately *physical*
        completeness notion, independent of the schema-format migration
        that :func:`find_data_gaps` tracks separately — an entry can be
        ``fully_defined`` whether it stores legs as the legacy ``legs[]``
        or the v3 ``trajectory.segments``.

        Dispatches by :attr:`cycler_class` (spec §16.7.5):

        * **single-ellipse**: ``orbit_elements.a_au`` + ``e`` present,
          V∞ at every encounter, at least one leg with ``tof_days``.
        * **multi-arc**: ``invariants`` dict present and non-empty
          (at least one non-null value), plus V∞ + tof_days guard
          (same leg/encounter check as single-ellipse).
        * **non-keplerian**: ``cr3bp`` identity triple
          (``jacobi_constant`` / ``period_nd`` / ``stability_index``)
          all non-null.  V∞ / legs guard intentionally *not* applied —
          Kepler elements are structurally inapplicable.
        """
        if self.data_gaps:
            return False
        if self.cycler_class == "non-keplerian":
            cr = self.cr3bp or {}
            return all(
                cr.get(k) is not None for k in ("jacobi_constant", "period_nd", "stability_index")
            )
        if self.cycler_class == "multi-arc":
            core_ok = bool(self.invariants) and any(
                v is not None for v in (self.invariants or {}).values()
            )
        else:  # single-ellipse (default)
            core_ok = self.orbit_elements_a_au is not None and self.orbit_elements_e is not None
        if not core_ok:
            return False
        if not self.vinf_kms_at_encounters or any(
            v is None for _, v in self.vinf_kms_at_encounters
        ):
            return False
        return bool(self.legs_tof_days) and all(t is not None for t in self.legs_tof_days)


@dataclass(frozen=True)
class MatchResult:
    """Tagged result of :func:`match` per spec §16.3."""

    outcome: Literal["known", "probable-match-NEEDS-HUMAN", "novel"]
    entry: CatalogueEntry | None
    distance: float | None


# ---------------------------------------------------------------------------
# Private helpers — sequence canonicalisation + binning + sha1
# ---------------------------------------------------------------------------


def _lex_min_rotation(sequence_str: str) -> str:
    """Return the lexicographically-minimal cyclic rotation of a body sequence.

    Operates on dash-joined body code lists. ``"E-M-E"`` → ``"E-E-M"``;
    ``"E-M"`` → ``"E-M"`` (already minimal); ``"E"`` → ``"E"``.

    Algorithm: simple O(n²) comparison of every rotation. For
    ``L_max = 8`` this is ~64 string comparisons per call, negligible
    in the M7 budget. Booth's O(n) algorithm is a future-work swap-in
    if profiling ever shows this as a hot path; M8's VEM at
    ``L_max <= 12`` is still well below the breakpoint.
    """
    parts = sequence_str.split("-")
    n = len(parts)
    if n <= 1:
        return sequence_str
    best = parts
    for i in range(1, n):
        rotated = parts[i:] + parts[:i]
        if rotated < best:
            best = rotated
    return "-".join(best)


def _dedupe_closing_body(bodies: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    """Drop the trailing body if it equals the first (open-sequence convention).

    M5's idealised cyclers store the closed loop as an open sequence
    that ends at the starting body — e.g. ``("E", "M", "E")`` for a
    2-syn E-M-E cell. The catalogue's ``sequence_canonical`` stores
    the same identity as ``"E-M"``. Drop the trailing repeat so the
    two representations canonicalise to the same string.
    """
    tup = tuple(bodies)
    if len(tup) >= 2 and tup[0] == tup[-1]:
        return tup[:-1]
    return tup


def _bin_vinf(vinf_kms: float) -> float:
    """Round V∞ to the nearest :data:`VINF_BIN_KMS` bin centre.

    Python's ``round`` uses banker's rounding; spec §16.2 does not pin
    a rounding rule and banker's-rounding ties never occur in
    practice on real V∞ values (0.025 km/s is below the precision of
    every catalogue source).
    """
    return round(vinf_kms / VINF_BIN_KMS) * VINF_BIN_KMS


def _bin_a_au(a_au: float) -> float:
    """Round semi-major axis (AU) to the nearest :data:`A_BIN_AU` bin centre."""
    return round(a_au / A_BIN_AU) * A_BIN_AU


def _bin_e(e: float) -> float:
    """Round eccentricity to the nearest :data:`E_BIN` bin centre."""
    return round(e / E_BIN) * E_BIN


def _canonical_json(d: dict[str, Any]) -> str:
    """Sorted-key, no-whitespace, ASCII-only JSON serialisation."""
    return json.dumps(d, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _signature_hash(input_dict: dict[str, Any]) -> str:
    """Compute the ``"sha1:..."`` digest of the canonical JSON of ``input_dict``."""
    blob = _canonical_json(input_dict).encode("utf-8")
    return "sha1:" + hashlib.sha1(blob).hexdigest()


def _assemble_signature(
    *,
    raw_bodies: tuple[str, ...] | list[str],
    raw_sequence: str,
    sense: SenseT,
    period_k: int,
    period_years: float,
    vinf_multiset_raw: list[tuple[str, float]],
    leg_elements_raw: list[tuple[float, float]],
    model_assumption: str,
) -> CanonicalSignature:
    """Single canonicalisation core for spec §16.2 signatures.

    THE one and only place that bins, sorts, assembles the hash-input
    dict, and computes the digest. Both signature entry points —
    :func:`canonical_signature` (recomputed from a model
    :class:`Cycler`) and :func:`_signature_from_yaml_fields` (published
    catalogue values) — extract their raw ``(bodies, sequence, V∞
    multiset, (a, e) multiset)`` from their respective sources and then
    delegate here, so a single identity can never hash two ways.

    All inputs are **raw** (un-binned): this function applies
    :func:`_dedupe_closing_body`, :func:`_lex_min_rotation`, and the
    ``_bin_*`` quantisers. ``period_years`` and ``model_assumption`` are
    carried on the returned :class:`CanonicalSignature` but, per the
    user-resolved decision #2 and plan §5 risk #15, are NOT in the hash.
    """
    deduped = _dedupe_closing_body(raw_bodies)
    bodies_sorted = tuple(sorted(set(deduped)))
    sequence_min = _lex_min_rotation(raw_sequence)
    vinf_binned = sorted((body, _bin_vinf(float(v))) for body, v in vinf_multiset_raw)
    # Dedupe per the per-cycler-signature convention: cyclers whose
    # multiple "legs" trace arcs of the *same* heliocentric ellipse
    # (e.g. the Aldrin family) must land on the same multiset as a
    # catalogue row that broadcasts one orbit_elements block over
    # several legs. A sorted set collapses both representations.
    leg_binned = sorted({(_bin_a_au(float(a)), _bin_e(float(e))) for a, e in leg_elements_raw})
    hash_input: dict[str, Any] = {
        "bodies": list(bodies_sorted),
        "sequence_canonical": sequence_min,
        "sense": sense,
        "period_k": int(period_k),
        "vinf_multiset_binned": [list(t) for t in vinf_binned],
        "leg_elements_multiset_binned": [list(t) for t in leg_binned],
    }
    return CanonicalSignature(
        bodies=bodies_sorted,
        sequence_canonical=sequence_min,
        sense=sense,
        period_k=int(period_k),
        period_years=float(period_years),
        vinf_multiset_binned=tuple(vinf_binned),
        leg_elements_multiset_binned=tuple(leg_binned),
        model_assumption=model_assumption,
        hash=_signature_hash(hash_input),
    )


# ---------------------------------------------------------------------------
# Public signature computation
# ---------------------------------------------------------------------------


def canonical_signature(
    cycler: Cycler,
    *,
    model_assumption: str = "circular-coplanar",
    period_k: int | None = None,
    period_years: float | None = None,
) -> CanonicalSignature:
    """Compute the spec §16.2 canonical signature of an idealised cycler.

    Parameters
    ----------
    cycler:
        The M3-shipped :class:`Cycler` instance. Must have populated
        ``bodies``, ``encounters``, and ``legs``.
    model_assumption:
        Pool-filter tag — default ``"circular-coplanar"`` (the M5
        idealised regime). Not included in the hash; consumed by
        :meth:`Catalog.filter` at match time.
    period_k:
        Synodic-multiple integer for the cycler's period. M5's
        :class:`OptimisationResult` carries the cell ``.period_k``; if
        the caller does not pass this in, we default to ``1`` to keep
        the call site terse for tests that only need rotation /
        binning invariance.
    period_years:
        Period in Julian years for display. Not in the hash; consumers
        format this for catalogue display.

    Returns
    -------
    CanonicalSignature
        Frozen, with bitwise-stable :attr:`CanonicalSignature.hash`.
    """
    deduped = _dedupe_closing_body(cycler.bodies)

    # Raw (un-binned) V∞ multiset, one entry per encounter. Open-sequence
    # boundary encounters set vinf_in == vinf_out by construction (M3
    # convention); use vinf_in uniformly so the sum is rotation-stable.
    # On open sequences (e.g. ("E", "M", "E")) the closing body's
    # encounter duplicates the opening body's |V∞| — spec §16.2's
    # "multiset" is per-encounter, NOT per unique body.
    vinf_multiset_raw: list[tuple[str, float]] = [
        (enc.body, float(np.linalg.norm(enc.vinf_in))) for enc in cycler.encounters
    ]
    # Raw per-leg heliocentric (a, e), recomputed from model state. The
    # core dedupes/bins so the Aldrin family (multiple arcs of one
    # ellipse) lands on the same multiset as the catalogue's single
    # orbit_elements block (see _assemble_signature).
    leg_elements_raw: list[tuple[float, float]] = [
        tuple(orbit_elements_au(enc.r, leg.v_depart, MU_SUN_KM3_S2))  # type: ignore[misc]
        for leg, enc in zip(cycler.legs, cycler.encounters[:-1], strict=True)
    ]

    return _assemble_signature(
        raw_bodies=cycler.bodies,
        raw_sequence="-".join(deduped),
        sense=cycler.sense,
        period_k=1 if period_k is None else period_k,
        period_years=0.0 if period_years is None else period_years,
        vinf_multiset_raw=vinf_multiset_raw,
        leg_elements_raw=leg_elements_raw,
        model_assumption=model_assumption,
    )


def _signature_from_yaml_fields(
    bodies: tuple[str, ...],
    sequence_canonical: str,
    sense: SenseT,
    period_k: int,
    period_years: float,
    vinf_multiset_raw: list[tuple[str, float]],
    leg_elements_raw: list[tuple[float, float]],
    model_assumption: str,
) -> CanonicalSignature:
    """Compute a :class:`CanonicalSignature` directly from catalogue YAML fields.

    Catalogue entries store their published ``a, e, V∞`` values
    directly; this is the **source-data** entry point. It extracts the
    raw published multisets and delegates to :func:`_assemble_signature`
    — the single canonicalisation core shared with
    :func:`canonical_signature` — so the published-value path and the
    model-recomputed path can never hash the same identity two ways.
    """
    return _assemble_signature(
        raw_bodies=bodies,
        raw_sequence=sequence_canonical,
        sense=sense,
        period_k=period_k,
        period_years=period_years,
        vinf_multiset_raw=vinf_multiset_raw,
        leg_elements_raw=leg_elements_raw,
        model_assumption=model_assumption,
    )


# ---------------------------------------------------------------------------
# Catalogue loader
# ---------------------------------------------------------------------------


def _segments_as_legs(row: dict[str, Any]) -> list[dict[str, Any]]:
    """Return the leg-like list for signature / ToF / rev extraction.

    Schema v3 (spec §16.6.2) moves the per-leg conic arcs into
    ``trajectory.segments`` (OCM ``TRAJ``). When present that is the
    authoritative leg source; otherwise the legacy top-level ``legs[]``
    is used. Segments reuse the ``tof_days`` / ``n_revs`` keys, so the
    downstream loops are agnostic to which form supplied the list.
    """
    segments = (row.get("trajectory") or {}).get("segments")
    if segments:
        return list(segments)
    return list(row.get("legs") or [])


def _extract_signature_inputs(
    entry_dict: dict[str, Any],
) -> tuple[list[tuple[str, float]], list[tuple[float, float]]] | None:
    """Pull the V∞ multiset + leg ``(a, e)`` multiset out of a YAML row.

    Returns ``None`` only when the V∞ multiset is unavailable —
    without V∞ values the matcher has nothing to compare against
    (family-seed and citation-only rows).

    Entries lacking ``orbit_elements.a_au`` / ``e`` (most Russell /
    McConaghy table rows) still produce a signature with an EMPTY
    leg-elements multiset; the matcher pool filter and the V∞
    multiset still partition the identity space adequately for these
    rows. The empty-multiset signatures are consumable; what they
    cannot do is differentiate two rows whose only difference is
    ``(a, e)``.
    """
    vinfs_raw = entry_dict.get("vinf_kms_at_encounters") or []
    if not vinfs_raw:
        return None
    vinf_pairs: list[tuple[str, float]] = []
    for entry in vinfs_raw:
        body = entry.get("body")
        v = entry.get("vinf_kms")
        if body is None or v is None:
            return None
        vinf_pairs.append((str(body), float(v)))

    orbit = entry_dict.get("orbit_elements") or {}
    a_au = orbit.get("a_au")
    e = orbit.get("e")
    legs_raw = _segments_as_legs(entry_dict)
    if a_au is None or e is None or not legs_raw:
        # Signature is still computable from V∞ + sequence + (period_k,
        # bodies, sense). Leave the leg-elements multiset empty.
        return vinf_pairs, []
    # Catalogue entries currently tabulate a single (a, e) on the
    # orbital-elements block; broadcast that over every leg, then the
    # signature canonicalisation dedupes to the unique set. For the
    # Aldrin family this collapses two legs of the same ellipse to one
    # multiset entry.
    leg_elements: list[tuple[float, float]] = [(float(a_au), float(e)) for _ in legs_raw]
    return vinf_pairs, leg_elements


def _entry_from_yaml(row: dict[str, Any]) -> CatalogueEntry:
    """Project one YAML row into a :class:`CatalogueEntry`.

    Computes :attr:`signature` / :attr:`signature_hash` when the row
    has the structural fields :func:`_extract_signature_inputs` needs;
    otherwise leaves them ``None`` so the entry never participates in
    matching.
    """
    bodies = tuple(row.get("bodies") or [])
    sequence_canonical = str(row.get("sequence_canonical") or "")
    sense_raw = row.get("sense") or "n/a"
    sense: SenseT
    if sense_raw == "outbound":
        sense = "outbound"
    elif sense_raw == "inbound":
        sense = "inbound"
    else:
        sense = "n/a"
    period = row.get("period") or {}
    period_k_raw = period.get("k")
    period_k = int(period_k_raw) if period_k_raw is not None else None
    period_years_raw = period.get("years")
    period_years = float(period_years_raw) if period_years_raw is not None else None
    period_pair = period.get("pair")
    model_assumption = str(row.get("model_assumption") or "circular-coplanar")

    vinfs_raw = row.get("vinf_kms_at_encounters") or []
    vinf_pairs: list[tuple[str, float | None]] = []
    for v_entry in vinfs_raw:
        body = v_entry.get("body")
        v = v_entry.get("vinf_kms")
        vinf_pairs.append(
            (str(body) if body is not None else "?", float(v) if v is not None else None)
        )

    legs_raw = _segments_as_legs(row)
    leg_tofs: list[float | None] = []
    leg_revs: list[int] = []
    for leg in legs_raw:
        tof = leg.get("tof_days")
        leg_tofs.append(float(tof) if tof is not None else None)
        leg_revs.append(int(leg.get("n_revs") or 0))

    oe = row.get("orbit_elements") or {}
    a_au_raw = oe.get("a_au")
    e_raw = oe.get("e")
    i_raw = oe.get("inclination_deg")

    # Schema v4 fields (spec §16.7, 2026-06-03)
    cycler_class = str(row.get("cycler_class") or "single-ellipse")
    orbit_elements_reference_frame = str(oe.get("reference_frame") or "heliocentric-inertial")
    orbit_elements_center = str(oe.get("center") or "Sun")
    invariants_raw = row.get("invariants")
    invariants: dict[str, Any] | None = (
        dict(invariants_raw) if isinstance(invariants_raw, dict) else None
    )
    cr3bp_raw = oe.get("cr3bp")
    cr3bp: dict[str, Any] | None = dict(cr3bp_raw) if isinstance(cr3bp_raw, dict) else None
    basis_raw = period.get("basis")
    period_basis: tuple[dict[str, Any], ...] | None = (
        tuple(dict(b) for b in basis_raw) if basis_raw else None
    )
    # Schema v4.1 fields (spec §16.7.7, 2026-06-03) — arc-type descriptors
    free_return_arcs_raw = row.get("free_return_arcs")
    free_return_arcs: tuple[dict[str, Any], ...] | None = None
    if isinstance(free_return_arcs_raw, list):
        free_return_arcs = tuple(dict(a) for a in free_return_arcs_raw)

    sig: CanonicalSignature | None = None
    sig_hash: str | None = None
    inputs = _extract_signature_inputs(row)
    if (
        inputs is not None
        and period_k is not None
        and period_years is not None
        and len(bodies) >= 1
        and sequence_canonical
    ):
        vinf_pairs_input, leg_inputs = inputs
        sig = _signature_from_yaml_fields(
            bodies=bodies,
            sequence_canonical=sequence_canonical,
            sense=sense,
            period_k=period_k,
            period_years=period_years,
            vinf_multiset_raw=vinf_pairs_input,
            leg_elements_raw=leg_inputs,
            model_assumption=model_assumption,
        )
        sig_hash = sig.hash

    priority_date = row.get("priority_date")
    return CatalogueEntry(
        id=str(row.get("id") or ""),
        name=str(row.get("name") or row.get("id") or ""),
        source=str(row.get("source") or "literature"),
        trajectory_regime=str(row.get("trajectory_regime") or "ballistic"),
        model_assumption=model_assumption,
        bodies=bodies,
        sequence_canonical=sequence_canonical,
        sense=sense,
        period_pair=str(period_pair) if period_pair is not None else None,
        period_k=period_k,
        period_years=period_years,
        vinf_kms_at_encounters=tuple(vinf_pairs),
        legs_tof_days=tuple(leg_tofs),
        legs_n_revs=tuple(leg_revs),
        orbit_elements_a_au=float(a_au_raw) if a_au_raw is not None else None,
        orbit_elements_e=float(e_raw) if e_raw is not None else None,
        orbit_elements_i_deg=float(i_raw) if i_raw is not None else None,
        primary=str(row.get("primary") or "Sun"),
        first_published=row.get("first_published"),
        priority_date=str(priority_date) if priority_date is not None else None,
        our_status=row.get("our_status"),
        signature=sig,
        signature_hash=sig_hash,
        validation=dict(row.get("validation") or {}),
        discovery=dict(row.get("discovery") or {}),
        discovery_run=row.get("discovery_run"),
        raw=row,
        family=row.get("family"),
        data_gaps=tuple(row.get("data_gaps") or ()),
        cycler_class=cycler_class,
        orbit_elements_reference_frame=orbit_elements_reference_frame,
        orbit_elements_center=orbit_elements_center,
        invariants=invariants,
        cr3bp=cr3bp,
        period_basis=period_basis,
        free_return_arcs=free_return_arcs,
    )


@dataclass(frozen=True)
class Catalog:
    """Loaded, indexed projection of ``data/catalogue.yaml`` (spec §16.3 pool).

    Attributes
    ----------
    entries:
        Tuple of :class:`CatalogueEntry` in file-declaration order.
    by_id:
        ``id`` → :class:`CatalogueEntry` index.
    by_hash:
        :attr:`CanonicalSignature.hash` → :class:`CatalogueEntry`
        index. Family-seed / citation-only entries (whose
        :attr:`CatalogueEntry.signature_hash` is ``None``) are
        excluded.
    """

    entries: tuple[CatalogueEntry, ...]
    by_id: dict[str, CatalogueEntry]
    by_hash: dict[str, CatalogueEntry]

    def filter(
        self,
        *,
        bodies: tuple[str, ...] | None = None,
        k: int | None = None,
        model_assumption: str | None = None,
    ) -> tuple[CatalogueEntry, ...]:
        """Spec §16.3 matcher pool pre-filter.

        Returns the entries that share the supplied ``bodies`` (set
        equality after dedup), ``k``, and ``model_assumption``.
        ``None`` arguments skip the corresponding constraint.
        """
        out: list[CatalogueEntry] = []
        wanted_bodies = set(_dedupe_closing_body(bodies)) if bodies is not None else None
        for entry in self.entries:
            entry_bodies_set = set(_dedupe_closing_body(entry.bodies))
            if wanted_bodies is not None and entry_bodies_set != wanted_bodies:
                continue
            if k is not None and entry.period_k != k:
                continue
            if model_assumption is not None and entry.model_assumption != model_assumption:
                continue
            out.append(entry)
        return tuple(out)


def load_catalog(path: Path | str = CATALOGUE_PATH) -> Catalog:
    """Read ``data/catalogue.yaml`` and build a :class:`Catalog`.

    Computes the canonical signature for every entry whose YAML row
    carries the required structural fields. Family-seed and
    citation-only entries are retained in :attr:`Catalog.entries`
    with :attr:`CatalogueEntry.signature` set to ``None``.

    Parameters
    ----------
    path:
        Override for the on-disk catalogue path. Defaults to
        :data:`CATALOGUE_PATH`. Tests can point at an alternative
        catalogue without monkey-patching.

    Returns
    -------
    Catalog
        Frozen, indexed.
    """
    raw = yaml.safe_load(Path(path).read_text())
    entries: list[CatalogueEntry] = []
    by_id: dict[str, CatalogueEntry] = {}
    by_hash: dict[str, CatalogueEntry] = {}
    for row in raw:
        entry = _entry_from_yaml(row)
        entries.append(entry)
        if entry.id:
            by_id[entry.id] = entry
        if entry.signature_hash is not None and entry.signature_hash not in by_hash:
            # First-wins on hash collisions: Russell-table rows whose
            # YAML lacks ``orbit_elements.a_au`` / ``e`` share the
            # same ``(bodies, k, V∞ multiset)`` signature inputs and
            # therefore collide on hash. The by_hash index keeps the
            # earliest-declared entry; the others remain reachable
            # via :attr:`by_id` and are returned to the matcher as
            # "known" against the first-declared cousin. This is a
            # catalogue-data-completeness limitation, not a matcher
            # bug — adding ``a, e`` to those YAML rows breaks the
            # collisions automatically.
            by_hash[entry.signature_hash] = entry
    return Catalog(
        entries=tuple(entries),
        by_id=by_id,
        by_hash=by_hash,
    )


def find_data_gaps(catalog: Catalog) -> list[dict[str, Any]]:
    """Enumerate every known-unknown across the catalogue (spec §16.6.4).

    Surfaces two gap kinds, each returned as a dict carrying the gap
    fields plus an ``entry_id`` key:

    - **Explicit** ``data_gaps[]`` entries declared on a row (with
      ``path`` / ``kind`` / ``note`` / ``source_hint`` / ``todo_ref``) —
      the machine-readable known-unknowns the lazy backfill consumes.
    - **Structural** migration gaps: an entry still carrying a legacy
      top-level ``legs[]`` rather than ``trajectory.segments`` (spec
      §16.6.5). An un-migrated entry is itself a tracked gap.

    The result is a query over the catalogue, replacing any manual
    gap audit — this is the sweep behind "identify gaps and lazy-fill".
    """
    gaps: list[dict[str, Any]] = []
    for entry in catalog.entries:
        for gap in entry.data_gaps:
            rec = dict(gap)
            rec["entry_id"] = entry.id
            gaps.append(rec)
        has_segments = bool((entry.raw.get("trajectory") or {}).get("segments"))
        if entry.raw.get("legs") and not has_segments:
            gaps.append(
                {
                    "entry_id": entry.id,
                    "path": "trajectory.segments",
                    "kind": "derive",
                    "note": "entry still on legacy legs[]; not migrated to trajectory{}",
                    "source_hint": None,
                    "todo_ref": "#65-backfill",
                }
            )
    return gaps


# ---------------------------------------------------------------------------
# Distance + match
# ---------------------------------------------------------------------------


def _multiset_l1(
    a: tuple[tuple[Any, ...], ...],
    b: tuple[tuple[Any, ...], ...],
    weights: tuple[float, ...],
) -> float:
    """L1 distance between two multisets of fixed-arity tuples.

    Pads the shorter multiset by treating missing elements as
    "infinity away" (contributes the heaviest element's full weight
    per remaining slot). For our use the two multisets always have
    the same cardinality on the pool-filtered match path (the pool
    filter shares ``(bodies, k, model_assumption)``), so the padding
    case is a defensive guard, not a hot path.
    """
    # Tuples are already sorted lexicographically by construction
    # (canonical_signature sorts them). Pair them positionally.
    n = max(len(a), len(b))
    total = 0.0
    for i in range(n):
        ta = a[i] if i < len(a) else None
        tb = b[i] if i < len(b) else None
        if ta is None or tb is None:
            # Missing element — count one full bin per position per
            # numeric slot.
            total += sum(weights)
            continue
        for j, w in enumerate(weights):
            # Slot j may be a string (body code) for the V∞ multiset
            # position 0; weight it 0 if equal, else one full bin.
            ja = ta[j]
            jb = tb[j]
            if isinstance(ja, str) or isinstance(jb, str):
                total += 0.0 if ja == jb else 1.0
            else:
                total += abs(float(ja) - float(jb)) / w
    return total


def signature_distance(
    sig_a: CanonicalSignature,
    sig_b: CanonicalSignature,
) -> float:
    """Weighted-L1 distance between two signatures per M7 plan §4.4.

    Distance accumulates:

    * ``|Δperiod_years| / 0.5`` (informational; ``period_years`` is
      not in the hash but the distance treats it as a soft
      tie-breaker).
    * ``Σ |Δvinf_kms| / 0.05`` over the V∞ multiset (paired by
      position after sorting).
    * ``Σ (|Δa_au| / 0.01 + |Δe| / 0.01)`` over the leg-elements
      multiset.

    A distance of ``0.0`` indicates per-position bin equality on every
    field; ``< TAU_NEAR`` indicates probable-match; ``>= TAU_NEAR``
    indicates novel.
    """
    period_term = abs(sig_a.period_years - sig_b.period_years) / 0.5
    vinf_term = _multiset_l1(
        sig_a.vinf_multiset_binned,
        sig_b.vinf_multiset_binned,
        weights=(0.0, VINF_BIN_KMS),
    )
    leg_term = _multiset_l1(
        sig_a.leg_elements_multiset_binned,
        sig_b.leg_elements_multiset_binned,
        weights=(A_BIN_AU, E_BIN),
    )
    return period_term + vinf_term + leg_term


def match(candidate: CanonicalSignature, catalog: Catalog) -> MatchResult:
    """Spec §16.3 three-state matcher: ``known`` / ``probable`` / ``novel``.

    Pool filter by ``(model_assumption, bodies, k)`` per spec §12.2 +
    §16.3. Exact hits are returned at distance ``0.0``; near hits
    below :data:`TAU_NEAR` are routed to ``probable-match-NEEDS-HUMAN``;
    everything else is ``novel``.
    """
    pool = catalog.filter(
        bodies=candidate.bodies,
        k=candidate.period_k,
        model_assumption=candidate.model_assumption,
    )
    # 1. Exact match (within binning) — hash lookup.
    by_hash = {e.signature_hash: e for e in pool if e.signature_hash is not None}
    if candidate.hash in by_hash:
        return MatchResult(outcome="known", entry=by_hash[candidate.hash], distance=0.0)
    # 2. Probable match — weighted L1 below TAU_NEAR.
    near: list[tuple[CatalogueEntry, float]] = []
    for entry in pool:
        if entry.signature is None:
            continue
        d = signature_distance(candidate, entry.signature)
        if d < TAU_NEAR:
            near.append((entry, d))
    if near:
        best_entry, best_distance = min(near, key=lambda x: x[1])
        return MatchResult(
            outcome="probable-match-NEEDS-HUMAN",
            entry=best_entry,
            distance=best_distance,
        )
    # 3. Novel.
    return MatchResult(outcome="novel", entry=None, distance=None)


# ---------------------------------------------------------------------------
# Frozen-immutability helper used by writeback module
# ---------------------------------------------------------------------------


def _replace_entry(entry: CatalogueEntry, **changes: Any) -> CatalogueEntry:
    """Return a new :class:`CatalogueEntry` with ``changes`` applied.

    Thin wrapper around :func:`dataclasses.replace` so writeback
    callers don't need to import ``dataclasses`` directly.
    """
    return replace(entry, **changes)


__all__ = [
    "A_BIN_AU",
    "CATALOGUE_PATH",
    "E_BIN",
    "TAU_NEAR",
    "VINF_BIN_KMS",
    "CanonicalSignature",
    "Catalog",
    "CatalogueEntry",
    "MatchResult",
    "canonical_signature",
    "find_data_gaps",
    "load_catalog",
    "match",
    "signature_distance",
]


def _main(argv: list[str] | None = None) -> int:
    """Minimal CLI: sweep the catalogue for known-unknowns (spec §16.6.4).

    ``python -m cyclerfinder.data.catalog gaps [--catalog PATH]`` prints
    the :func:`find_data_gaps` result as a JSON array.
    """
    import argparse

    parser = argparse.ArgumentParser(prog="cyclerfinder.data.catalog")
    sub = parser.add_subparsers(dest="cmd", required=True)
    gaps = sub.add_parser("gaps", help="enumerate data_gaps[] + migration gaps")
    gaps.add_argument("--catalog", default=str(CATALOGUE_PATH))
    args = parser.parse_args(argv)

    if args.cmd == "gaps":
        catalog = load_catalog(args.catalog)
        print(json.dumps(find_data_gaps(catalog), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
