"""Catalogue + identity subpackage (M7).

Hosts the spec §16 catalogue loader, canonical signature, identity
matcher, append-only ledger, and entry-mutation writeback helpers.

* :mod:`cyclerfinder.data.catalog` — frozen :class:`CatalogueEntry`,
  :class:`CanonicalSignature`, :func:`canonical_signature`,
  :class:`Catalog`, :func:`load_catalog`, :func:`match`,
  :class:`MatchResult`. Spec §16.1, §16.2, §16.3, §12.2.
* :mod:`cyclerfinder.data.ledger` — JSONL append-only persistence layer
  (spec §13.6, §13.8) for resumable / parallel-safe finder runs.
* :mod:`cyclerfinder.data.writeback` — :func:`apply_v0_v1_to_entry` /
  :func:`apply_v2_to_entry` / :func:`apply_v3_to_entry` /
  :func:`record_rediscovery` / :func:`register_discovery` /
  :func:`serialise_entry_yaml`. Spec §16.1 schema-v2 writeback.
* :mod:`cyclerfinder.data.discover` — ledger-backed runner that wraps
  :func:`cyclerfinder.search.optimize.find_cyclers` with per-cell
  persistence, signature computation, and auto-validation.

Importing this subpackage is cheap — the catalogue YAML is only read
when :func:`~cyclerfinder.data.catalog.load_catalog` is called, and the
M6a/M6b verifier dependencies are pulled in lazily by
:mod:`~cyclerfinder.data.discover`.

Plan: ``docs/phases/m7-catalogue-novelty-matching/plan.md``.
"""
