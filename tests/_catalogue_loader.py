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
* At least one leg with non-null ``tof_days``.

Family-seed and citation-only entries (most fields ``null``) are filtered
out by the non-null checks above.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "seed_cyclers.yaml"


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


def load_constructible_entries() -> list[CatalogueEntry]:
    """Return the filtered, sortable list of catalogue entries the v1
    optimiser can construct.

    Sorted ascending by ``id`` so pytest parametrisation produces a
    deterministic ordering across runs and across machines.
    """
    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    entries: list[CatalogueEntry] = []
    for row in raw:
        if row.get("trajectory_regime", "ballistic") != "ballistic":
            continue
        if row.get("primary", "Sun") != "Sun":
            continue
        bodies = row.get("bodies") or []
        if len(bodies) != 2:
            continue
        sequence_canonical = row.get("sequence_canonical") or ""
        if not _is_two_body_alternation(sequence_canonical, tuple(bodies)):
            continue
        period = row.get("period") or {}
        if period.get("k") is None or period.get("years") is None:
            continue
        vinfs_raw = row.get("vinf_kms_at_encounters") or []
        if len(vinfs_raw) != len(bodies):
            continue
        vinfs = [v.get("vinf_kms") for v in vinfs_raw]
        if any(v is None for v in vinfs):
            continue
        legs_raw = row.get("legs") or []
        leg_tofs = [leg.get("tof_days") for leg in legs_raw if leg.get("tof_days") is not None]
        if not leg_tofs:
            continue

        entries.append(
            CatalogueEntry(
                id=row["id"],
                name=row.get("name", row["id"]),
                bodies=tuple(bodies),
                sequence_canonical=sequence_canonical,
                period_k=int(period["k"]),
                period_years=float(period["years"]),
                vinf_targets_kms=tuple(float(v) for v in vinfs),
                leg_tofs_days=tuple(float(t) for t in leg_tofs),
            )
        )
    entries.sort(key=lambda e: e.id)
    return entries
