# #279 disposition — Finley et al. JSR Pluto moon-tour paper

**Date:** 2026-06-15
**Disposition:** **DISSOLVED.** The orbital-tour content is contained within the
already-characterized Stern-Tapley-Finley-Scherrer 2020 paper
(DOI `10.2514/1.A34658`, *J. Spacecraft & Rockets* 57(5):956). No separate
Finley-first-author publication appears to exist.

## Investigation summary
Multiple WebSearch + WebFetch attempts (this session):
- arc.aiaa.org / doi.org 403-paywalled on every direct DOI URL.
- AAS Astrodynamics Specialist Conference proceedings 2018-2020: no matching
  separately-attributed Finley paper found.
- Tiffany Finley's SwRI profile page: HTTP 404 (URL deprecated).
- Semantic Scholar / ResearchGate / ADS: 403 or empty result.
- arXiv preprint search: zero hits matching the title.
- SwRI press release, EurekAlert (577420), ScienceDaily (181024163636), PhysOrg
  (2018-10), Centauri Dreams (2018-10-25), SpaceFlight Insider — all 6 popular
  coverages describe Finley as the **software lead who designed the tour**, and
  describe the tour as content *inside* the Stern et al. 2020 paper, not a
  separate publication.

The closest candidate is one WebSearch result citing the title "An Orbital Tour
of Pluto and Its Moons" as a JSR paper by Finley T., Barth E., Howett C.,
Zangari A., Tapley M., Scherrer J., and Stern A. — but the DOI is not surfaced;
all direct DOI URLs 403; and the SwRI press release explicitly attributes the
"accepted by JSR" paper to Amanda Zangari (the 45-objects Kuiper Belt paper,
DOI `10.2514/1.A34329`), not to Finley.

## Risk assessment for the 12 Pluto gauntlet candidates
- The Stern et al. 2020 paper (already characterized in #269 + #270) describes
  a Charon-gravity-assist tour structurally distinct from the 12 candidates'
  patched-conic Lambert-leg multi-rev fingerprint.
- IF a separate Finley-first-author paper exists, it would be on the same AIAA
  venue with the same structural framing (the SwRI press release would have
  flagged a structurally different approach). The risk that such a hidden
  paper structurally matches the 12 candidates is **very low**.
- **#270's verdict stands**: 0 of 12 candidates literature-matched against the
  characterized Pluto-system spacecraft-trajectory literature.

## Recommendation for V5 sign-off
**Proceed.** The acquisition risk surfaced by #270 is mitigated by this
investigation. If a future deeper review (e.g. accessing AIAA paywall via
institutional account) surfaces a separate Finley paper, re-run the literature
check on the 12 candidates against its structural fingerprint. Until then,
treat the Stern et al. 2020 paper as the canonical anchor and proceed with
the V0-V5 gauntlet (#274).

## Source URL (for future direct access)
- AIAA: https://arc.aiaa.org/doi/10.2514/1.A34658 (paywalled)
- ADS: https://ui.adsabs.harvard.edu/abs/2020JSpRo..57..956S/abstract (paywalled
  full text)
- ResearchGate: publication 341139909 (paywalled / 403)

When institutional access is available, the Stern et al. 2020 paper (and any
companion Finley paper if it exists) should be filed at
`cyclers_pdf/papers/stern-tapley-finley-scherrer-2020-pluto-orbiter-kuiper-belt-explorer-jsr-doi-10.2514-1.A34658.pdf`.
