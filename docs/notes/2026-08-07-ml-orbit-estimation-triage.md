# Triage: two ML/orbit-estimation papers, both OUT-OF-SCOPE

Two papers uploaded by the user 2026-08-07, both landing in the same subfield — triaged together.

## Peng & Bai 2021, "Fusion of a machine learning approach and classical orbit predictions"

**Source:** H. Peng & X. Bai, *Acta Astronautica* 184 (2021) 222-240. DOI `10.1016/j.actaastro.2021.04.017`.
Rutgers University (Dept. of Mechanical and Aerospace Engineering). Filed at
`cyclers_pdf/papers/peng-bai-2021-fusion-machine-learning-classical-orbit-predictions-actaastro-184-222-doi-10.1016-j.actaastro.2021.04.017.pdf`
(md5 `f172571c7b719e4062fecba197b986a9`), text-layer, 19 pages.

**Verdict: OUT-OF-SCOPE — triaged, not digested.**

This is a Space Situational Awareness (SSA) / orbit-determination paper: it fuses a machine-learning
approach with an Extended Kalman Filter to improve the accuracy and uncertainty quantification of
**tracking/predicting the orbit of resident space objects** from observational data (collision
avoidance, catalog maintenance, sensor scheduling). This is orbit ESTIMATION/TRACKING, not orbit
DESIGN — a different subfield from this project's own scope (CR3BP/CRNBP/CCR4BP cycler and resonant-
orbit trajectory design and discovery). Matches this project's own existing "ML/surrogate/GNC
background" triage category (see `CORPUS_INDEX.md`'s dedicated table) rather than warranting a full
digest-and-citation-mining pass.

**Author note (checked, not assumed):** this is a DIFFERENT Hao Peng than the two already-corpused
"Peng-Xu" ERTBP papers (`peng-xu-2015-...`, `peng-xu-2017-...`) — different institution, different
co-author (Xiaoli Bai here vs. Xu on the ERTBP papers), different subfield entirely. A name match
alone would have been misleading; filed under `peng-bai-2021-...` to keep the two authors visually
distinct in the corpus.

## Caldas & Soares 2024, "Machine Learning in Orbit Estimation: a Survey"

**Source:** F. Caldas & C. Soares, arXiv:2207.08993v4 (2024). Filed at
`cyclers_pdf/papers/caldas-soares-2024-machine-learning-orbit-estimation-survey-arxiv-2207.08993.pdf`
(md5 `2c0d3a6e06abf0df665d5e577ada56a9`), text-layer (LaTeX/pdfTeX), 39 pages.

**Verdict: OUT-OF-SCOPE — triaged, not digested.** A survey of ML methods applied to orbit
ESTIMATION (the paper's own title states this precisely) — the same SSA/orbit-determination
subfield as the Peng & Bai paper above, not orbit design. Same triage category, same reasoning:
this project's own scope is trajectory design/discovery (CR3BP/CRNBP/CCR4BP cyclers and resonant
orbits), not tracking/estimating the orbits of existing objects from observational data. As a
*survey*, it could in principle be a useful pointer to other ML-for-orbits work if this project's
own scope ever expanded into estimation — noted here rather than silently discarded, but not
citation-mined this pass (mining a triaged-out-of-scope survey is not this project's own established
practice; citation-mining is reserved for on-topic digested papers).

**Registration:** `CORPUS_INDEX.md`'s "ML / surrogate / GNC background (triaged sweep)" table updated
in the same session, for both papers.
