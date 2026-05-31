"""Catalogue data for cyclerfinder.

Houses the seed catalogue of published cyclers (``seed_cyclers.yaml``) that
later milestones consume:

- **M3** golden tests will pull the Aldrin and McConaghy 2-synodic entries as
  fixtures asserted against the constructed cyclers.
- **M7** novelty matching will check finder hits against this catalogue (and
  against entries added by ongoing literature review) to tag re-derivations.

No loader is shipped yet; consumers (M3 / M7) will add YAML deserialisation
when they need it. The data file's contract is the markdown reference at
``docs/known-cyclers.md`` plus the ``signature_fields`` shape in
``docs/spec.md`` §16.1.
"""
