# Crossref bibliographic search returns the AIAA conference DOI instead of the journal DOI for 1960s-2010s papers

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: no

What happened: `query.bibliographic` for Gillespie & Ross 1967 (JSR), Sohn 1964 (JSR), Hollister & Prussing 1966, Patel 1998 and Anderson 2018 returned the 10.2514/6.19xx-NNNN conference record first. Several journal DOIs were never found.
Workaround: added the journal name, volume and year to the query and checked volume and pages in the returned record. Where the journal version was missing, recorded the conference DOI explicitly as "conference version".
Suggested fix: a small helper (e.g. scripts/crossref_check.py) that requires a container-title or volume match before printing CONFIRMED.

Disposition (main, 2026-10-05): Backlog: #964 (corpus tooling: Crossref helper requiring container/volume match).
