"""Catalogue loader sanity + signature-index coverage (M7).

* Loader returns every row in ``data/catalogue.yaml`` (currently
  ~219; the assertion is "at least 200" to allow non-disruptive
  catalogue growth without test churn — exact count is asserted by
  the row-count-anchored test below).
* Constructible entries (per :func:`tests._catalogue_loader.load_constructible_entries`)
  have non-``None`` :attr:`CatalogueEntry.signature_hash`.
* Family-seed / citation-only entries have ``signature_hash is None``
  by design — they live in the catalogue but cannot be matched
  against.
* Pool partitioning by ``model_assumption`` correctly excludes CR3BP
  entries from circular-coplanar filters.
"""

from __future__ import annotations

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import CATALOGUE_PATH, load_catalog
from tests._catalogue_loader import load_constructible_entries


def test_load_catalog_loads_all_yaml_rows() -> None:
    """Loader yields one :class:`CatalogueEntry` per YAML row."""
    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    cat = load_catalog()
    assert len(cat.entries) == len(raw)
    assert len(cat.entries) >= 200, (
        f"catalogue size shrank unexpectedly: {len(cat.entries)} rows. "
        f"If this is a deliberate prune, update the lower bound."
    )


def test_load_catalog_by_id_index_populated() -> None:
    """Every entry with an ``id`` appears in :attr:`Catalog.by_id`."""
    cat = load_catalog()
    for entry in cat.entries:
        if entry.id:
            assert cat.by_id[entry.id] is entry


def test_load_catalog_contains_aldrin_classic() -> None:
    """Spec §16.4 anchor: Aldrin classic outbound is row 1."""
    cat = load_catalog()
    assert "aldrin-classic-em-k1-outbound" in cat.by_id
    aldrin = cat.by_id["aldrin-classic-em-k1-outbound"]
    assert aldrin.priority_date == "1985-10-28"
    assert aldrin.bodies == ("E", "M")
    assert aldrin.period_k == 1
    assert aldrin.sense == "outbound"


def test_load_catalog_contains_s1l1() -> None:
    """M5 binding-gate target: the 2-syn S1L1 entry is loaded."""
    cat = load_catalog()
    assert "s1l1-2syn-em-cpom" in cat.by_id
    entry = cat.by_id["s1l1-2syn-em-cpom"]
    assert entry.period_k == 2
    assert entry.bodies == ("E", "M")


def test_constructible_entries_have_signature_hash() -> None:
    """Every entry returned by the M5 constructibility loader has a
    populated signature_hash."""
    cat = load_catalog()
    constructible = load_constructible_entries()
    missing: list[str] = []
    for c in constructible:
        entry = cat.by_id.get(c.id)
        if entry is None:
            missing.append(f"{c.id} (not in catalogue)")
        elif entry.signature_hash is None:
            missing.append(f"{c.id} (signature_hash is None)")
    assert not missing, f"constructible entries missing signature_hash: {missing}"


def test_family_seed_entries_have_null_signature() -> None:
    """Family-seed / citation-only entries are loaded but never
    indexed by hash."""
    cat = load_catalog()
    family_seeds = {
        "jones-2017-vem-triple-family",
        "wittal-2022-em-cycler-family",
    }
    for fid in family_seeds:
        if fid in cat.by_id:
            entry = cat.by_id[fid]
            assert entry.signature_hash is None, (
                f"family-seed {fid} unexpectedly has signature_hash {entry.signature_hash!r}"
            )


def test_catalog_filter_by_bodies() -> None:
    """``filter(bodies=("E","M"))`` returns only Earth-Mars entries."""
    cat = load_catalog()
    em_only = cat.filter(bodies=("E", "M"))
    assert len(em_only) > 0
    for entry in em_only:
        # Body set after dedup (closing-body convention) must equal {"E","M"}.
        body_set = set(entry.bodies)
        if entry.bodies and entry.bodies[0] == entry.bodies[-1]:
            body_set = set(entry.bodies[:-1])
        assert body_set == {"E", "M"}, (
            f"filter bodies=(E,M) leaked non-E/M entry: {entry.id} bodies={entry.bodies}"
        )


def test_catalog_filter_by_k() -> None:
    """``filter(k=1)`` returns only ``period_k == 1`` entries."""
    cat = load_catalog()
    k1 = cat.filter(k=1)
    for entry in k1:
        assert entry.period_k == 1, f"filter k=1 leaked entry {entry.id} with k={entry.period_k}"


def test_catalog_filter_by_model_assumption_partitions_cr3bp() -> None:
    """**Spec §12.2 partitioning binding gate.**

    CR3BP entries (Arenstorf, etc.) live in a different pool from
    circular-coplanar; filtering by one excludes the other.
    """
    cat = load_catalog()
    cr3bp = cat.filter(model_assumption="cr3bp")
    cpom = cat.filter(model_assumption="circular-coplanar")
    cr3bp_ids = {e.id for e in cr3bp}
    cpom_ids = {e.id for e in cpom}
    # Aldrin classic is circular-coplanar; should NOT appear in CR3BP filter.
    assert "aldrin-classic-em-k1-outbound" in cpom_ids
    assert "aldrin-classic-em-k1-outbound" not in cr3bp_ids
    # No overlap.
    assert not (cr3bp_ids & cpom_ids), f"model_assumption pools overlap: {cr3bp_ids & cpom_ids}"


def test_load_catalog_deterministic() -> None:
    """Two consecutive loads produce equal :attr:`Catalog.by_hash`."""
    a = load_catalog()
    b = load_catalog()
    assert a.by_hash == b.by_hash
    assert len(a.entries) == len(b.entries)
