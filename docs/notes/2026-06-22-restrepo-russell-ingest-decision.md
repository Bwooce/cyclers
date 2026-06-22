# Decision Note: Restrepo-Russell JPL 3BP Catalogue Ingest
Date: 2026-06-22
Context: Task #377 (Restrepo-Russell 2018 Phase 2/3 Ingest Decision)

## Decision: Defer Wholesale Ingest

We have decided to **DEFER** the wholesale ingestion of the Restrepo-Russell ~3 million periodic orbit database. Instead, we will selectively ingest the *connecting resonance* files (`H_R p:q - LL1/LL2/DPO/QDRO`) only when the discovery daemon actually pushes a planar CR3BP candidate that requires it as an establishment seed or a novelty cross-check.

## Justification for Ignoring Bound Libration Families and Pure Resonances

The vast majority of the 3 million orbits in the database are not cycler material in our project's `orbit_class` taxonomy. Specifically:

*   **Bound Libration Families (LL1, LL2, DRO, DPO, Hg, Hb, etc.):** These orbits remain in the vicinity of the secondary and never leave its sphere of influence. They loiter rather than cycle, making them fundamentally incompatible with our cycler-search focus.
*   **Pure Resonances (R_{p:q}):** Most simple resonances orbit only the primary without ever having a secondary encounter. While some inner-resonance cases with close approaches might be precursors, they are not cyclers.
*   **Catalogue Inflation and V2 Failure:** The project's sourced-cycler bar is high. Ingesting these non-cycler families wholesale would massively inflate the catalogue with rows that fail the V2 validation gate (the cycler-definition gate) by construction, diluting the quality of the discovery sweeps.

## Next Steps

1. The `CorpusAnchor` for the database has been added to `KNOWN_CORPUS` in `literature_check.py`. This ensures that any CR3BP planar-axisymmetric candidate hitting this structural footprint will be correctly flagged as a "published" rediscovery.
2. We will implement a focused ingest path (`scripts/ingest_restrepo_russell_database.py`) if and when the first discovery-daemon SILVER candidate in a Restrepo-Russell-covered system reaches a Jacobi constant inside the database's energy envelope.
