"""Load the seed catalogue, filter to entries the v1 optimiser can construct.

Used by ``tests/test_catalogue_rediscovery.py`` to parametrise the
rediscovery test over the catalogue. As the catalogue grows (more Russell
rows, more McConaghy SnLm entries, finder discoveries, etc.) the test
suite grows automatically.

Filter rules — an entry is *constructible by the v1 optimiser* iff:

* ``trajectory_regime == "ballistic"`` (absent → defaults to ``"ballistic"``).
* ``primary == "Sun"`` (absent → defaults to ``"Sun"``).
* ``len(bodies) == 2`` — M5's free-parameter dimension assumes a single
  synodic pair; 3-body / VEM optimisation is M8 work.
* ``sequence_canonical`` is a strict 2-body alternation across exactly two
  encounters (e.g. ``"E-M"`` or ``"M-E"``). The v1 optimiser builds the
  corresponding M4 cell with sequence ``(body_a, body_b, body_a)`` and
  ``period_k`` from the entry. Entries with longer canonical sequences
  (e.g. ``"E-E-M-M"`` for 2-synodic McConaghy / Russell rows) need M4's
  structural enumeration to find the matching cell and are excluded
  here — they're constructible *in principle*, but wiring them up
  requires extra structural inference outside the v1 loader's contract.
* All entries in ``vinf_kms_at_encounters`` have non-null ``vinf_kms``,
  and ``len(vinf_kms_at_encounters) == len(bodies)``.
* ``period.k`` and ``period.years`` are non-null.
* At least one leg with non-null ``tof_days`` (read from
  ``trajectory.segments`` when migrated, else the legacy ``legs[]``).

Family-seed and citation-only entries (most fields ``null``) are filtered
out by the non-null checks above.
"""

from __future__ import annotations

import dataclasses
import enum
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"


class ExclusionReason(enum.Enum):
    """Why a catalogue row is (or is not) admitted to the v1 rediscovery gauntlet.

    Exactly one reason applies to every catalogue row. ``CONSTRUCTIBLE``
    rows go into ``test_catalogue_entry_rediscovers``; every other reason
    documents a *categorised* exclusion so no entry vanishes silently as
    the catalogue grows (the coverage-audit invariant, task #55).

    The order of the non-constructible members mirrors the order of the
    filter checks in :func:`classify_row` — the *first* failing check
    assigns the reason.
    """

    CONSTRUCTIBLE = "constructible"
    """Passes every v1 filter: a heliocentric, ballistic/powered, 2-body
    alternation with non-null period, V∞, and at least one leg ToF. Eligible
    for V0→V1 promotion (the rediscovery test will run it, or skip it via
    a documented ``EXPECTED_SKIPS`` entry)."""

    NON_BALLISTIC = "non_ballistic_regime"
    """``trajectory_regime`` is neither ``ballistic`` nor ``powered`` (e.g.
    low-thrust Sims-Flanagan rows): not a Lambert cell the v1 optimiser builds."""

    NON_HELIOCENTRIC = "non_heliocentric"
    """``primary != "Sun"`` (lunar / Jovian moon tours). The v1 optimiser is
    heliocentric; these await planet-centric support (task #76)."""

    NOT_TWO_BODY = "not_two_body"
    """``len(bodies) != 2``: VEM triples and other multi-body itineraries need
    M4 structural enumeration over >1 synodic pair, outside the v1 contract."""

    MULTI_ENCOUNTER_SEQUENCE = "multi_encounter_sequence"
    """Two bodies, but ``sequence_canonical`` has more than two encounters
    (e.g. ``E-E-M-M`` for 2-synodic rows). Constructible *in principle* but
    needs multi-rev Lambert + structural inference the v1 loader does not do."""

    MISSING_PERIOD = "missing_period"
    """``period.k`` or ``period.years`` is null — typically a family-seed or
    citation-only row carrying no closed-orbit data."""

    MISSING_VINF = "missing_vinf"
    """``vinf_kms_at_encounters`` is absent, mis-sized, or has a null
    ``vinf_kms`` — no published target to rediscover against."""

    MISSING_LEG_TOFS = "missing_leg_tofs"
    """No leg carries a non-null ``tof_days`` (neither ``trajectory.segments``
    nor legacy ``legs[]``) — nothing to seed or compare leg geometry against."""


@dataclasses.dataclass(frozen=True)
class CatalogueEntry:
    """The minimal projection of a catalogue row needed for rediscovery.

    Attributes
    ----------
    id:
        Catalogue id (e.g. ``"aldrin-classic-em-k1-outbound"``).
    name:
        Human-readable name; defaults to ``id`` if absent on the YAML row.
    bodies:
        Two-tuple of body codes (e.g. ``("E", "M")``). The loader
        guarantees ``len(bodies) == 2``.
    sequence_canonical:
        Canonical cyclic sequence per spec §16.2 (e.g. ``"E-M"``). The
        loader guarantees this is a 2-body alternation matching
        ``bodies``.
    period_k:
        Period in synodic multiples (``period.k`` in the YAML).
    period_years:
        Period in Julian years (``period.years`` in the YAML).
    vinf_targets_kms:
        Tuple of published V∞ magnitudes at each encounter, in the same
        order as ``bodies``. The loader guarantees one entry per body
        and all non-null.
    leg_tofs_days:
        Tuple of non-null published leg times-of-flight in days, in leg
        order. Carried for diagnostics / future tof-comparison checks.
    """

    id: str
    name: str
    bodies: tuple[str, ...]
    sequence_canonical: str
    period_k: int
    period_years: float
    vinf_targets_kms: tuple[float, ...]
    leg_tofs_days: tuple[float, ...]


def _is_two_body_alternation(
    sequence_canonical: str,
    bodies: tuple[str, ...],
) -> bool:
    """Return True iff ``sequence_canonical`` is a 2-encounter alternation
    of the two ``bodies``.

    The v1 optimiser's cell template is
    ``Cell.sequence = (body_a, body_b, body_a)`` with ``period_k = K``;
    that matches catalogue rows whose canonical-cyclic form is two
    encounters (one of each body) — i.e. ``"E-M"`` or ``"M-E"`` for
    Earth-Mars rows.
    """
    if len(bodies) != 2:
        return False
    parts = [p for p in sequence_canonical.split("-") if p]
    if len(parts) != 2:
        return False
    return set(parts) == set(bodies)


def classify_row(row: dict[str, Any]) -> tuple[ExclusionReason, CatalogueEntry | None]:
    """Classify one catalogue row for the v1 rediscovery gauntlet.

    Single source of truth for both :func:`load_constructible_entries`
    (which keeps the ``CONSTRUCTIBLE`` rows) and the coverage-audit census
    (which tallies every reason). Returns the assigned
    :class:`ExclusionReason` and, when ``CONSTRUCTIBLE``, the parsed
    :class:`CatalogueEntry`; otherwise ``None``.

    The checks run in a fixed order and the *first* failure wins, so the
    reason is stable and the categories are mutually exclusive.
    """
    # Regime describes the ΔV character (ballistic vs gravity-assist-
    # maintained powered), NOT constructibility: both are buildable as a
    # Lambert cell here. Mirrors the m6b loader, which also admits both.
    if row.get("trajectory_regime", "ballistic") not in ("ballistic", "powered"):
        return ExclusionReason.NON_BALLISTIC, None
    if row.get("primary", "Sun") != "Sun":
        return ExclusionReason.NON_HELIOCENTRIC, None
    bodies = row.get("bodies") or []
    if len(bodies) != 2:
        return ExclusionReason.NOT_TWO_BODY, None
    sequence_canonical = row.get("sequence_canonical") or ""
    if not _is_two_body_alternation(sequence_canonical, tuple(bodies)):
        return ExclusionReason.MULTI_ENCOUNTER_SEQUENCE, None
    period = row.get("period") or {}
    if period.get("k") is None or period.get("years") is None:
        return ExclusionReason.MISSING_PERIOD, None
    vinfs_raw = row.get("vinf_kms_at_encounters") or []
    if len(vinfs_raw) != len(bodies):
        return ExclusionReason.MISSING_VINF, None
    vinfs = [v.get("vinf_kms") for v in vinfs_raw]
    if any(v is None for v in vinfs):
        return ExclusionReason.MISSING_VINF, None
    # Schema v3 (spec §16.6.2): per-leg arcs may live under
    # trajectory.segments (OCM TRAJ) once an entry is migrated off the
    # legacy flat legs[]. Mirror the package loader's _segments_as_legs
    # fallback so the harness reads both on-disk forms during the lazy
    # backfill (e.g. the migrated Aldrin classic out/in entries).
    legs_raw = (row.get("trajectory") or {}).get("segments") or row.get("legs") or []
    leg_tofs = [leg.get("tof_days") for leg in legs_raw if leg.get("tof_days") is not None]
    if not leg_tofs:
        return ExclusionReason.MISSING_LEG_TOFS, None

    entry = CatalogueEntry(
        id=row["id"],
        name=row.get("name", row["id"]),
        bodies=tuple(bodies),
        sequence_canonical=sequence_canonical,
        period_k=int(period["k"]),
        period_years=float(period["years"]),
        vinf_targets_kms=tuple(float(v) for v in vinfs),
        leg_tofs_days=tuple(float(t) for t in leg_tofs),
    )
    return ExclusionReason.CONSTRUCTIBLE, entry


def classify_catalogue() -> list[tuple[str, ExclusionReason]]:
    """Classify *every* catalogue row, returning ``(id, reason)`` sorted by id.

    Coverage-audit primitive (task #55): no row is dropped silently, so the
    census test can assert the union of all reasons covers the full catalogue.
    """
    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    out = [(row["id"], classify_row(row)[0]) for row in raw]
    out.sort(key=lambda pair: pair[0])
    return out


def load_constructible_entries() -> list[CatalogueEntry]:
    """Return the filtered, sortable list of catalogue entries the v1
    optimiser can construct.

    Sorted ascending by ``id`` so pytest parametrisation produces a
    deterministic ordering across runs and across machines.
    """
    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    entries: list[CatalogueEntry] = []
    for row in raw:
        reason, entry = classify_row(row)
        if reason is ExclusionReason.CONSTRUCTIBLE and entry is not None:
            entries.append(entry)
    entries.sort(key=lambda e: e.id)
    return entries
