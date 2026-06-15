"""Track-A genome modules: parametric orbit families with topological invariants.

The genome package collects compact, parametric descriptions of orbit families
where the family member is selected by a small set of integer / continuous
labels (the "genome"). Each module exposes a sourced IC table, a reproducer that
re-flies the published seed, and a topology classifier that confirms the
recovered orbit lies in the right family.

Members of this package are *opt-in*; they re-use the dynamical primitives in
:mod:`cyclerfinder.core` and the correctors in :mod:`cyclerfinder.search` and do
not modify them. The package name is reserved by the project spec; see
:doc:`docs/spec.md` Track A.
"""
