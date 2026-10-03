# Digest - Grabowski, Bellome & Felicetti (2026), "A diversity-based strategy for asteroid tour design"

**Digested:** 2026-10-03 (text-layer PDF; all 15 pages read). **Citation:** Jan Grabowski, Andrea Bellome, Leonard Felicetti
(Cranfield University), "A diversity-based strategy for asteroid tour design", Advances in Space Research 77 (2026)
7039-7053, **DOI 10.1016/j.asr.2026.01.075**; received 1 April 2025, revised 7 January 2026, accepted 21 January 2026, online
27 January 2026; open access, CC BY. Extends the authors' 2024 IEEE Aerospace Conference paper (stated in the text). Filed in
the private paper corpus as
`grabowski-bellome-felicetti-2026-diversity-based-strategy-asteroid-tour-design-asr-77-7039-doi-10.1016-j.asr.2026.01.075.pdf`
(md5 `701469f2bbc107bf478b832a4ad873b0`, 15 pages).

## What it is
A combinatorial-search paper for multi-asteroid rendezvous tours. Standard beam search (BS) ranks partial tours by delta-v and
keeps the best BW (beam width); the authors observe it returns many tours with the same bodies in different epochs. They add
a diversity measure borrowed from the Maximum Diversity Problem and use it to choose which partial tours survive each level.

## Diversity score and search strategy (Q1)
- Diversity distance between two equal-length tours A and B (Eq. 6): Div(A,B) = [(dv(A)+dv(B))/2] * [(L - nb_common)/L],
  with L the number of bodies and nb_common the number of bodies in common. Identical tours, or the same bodies in a different
  order, give zero.
- Diversity score of a tour (Eq. 7): J_div(A_i) = sum over k other than i of Div(A_i, A_k), divided by (m-1): its average distance to
  the other tours in the set. It is relative to the set, so "the diversity scores of the same sequence in two different set of
  solutions cannot be compared". It deliberately includes total delta-v, on the stated assumption that "high Dv tours are more
  likely to improve the diversity of the solution set".
- Strategies, all with the same tree (nodes are object-at-epoch, Lambert-arc branches): BS keeps the BW lowest-delta-v tours;
  Diversity Search (DS) first removes duplicates (same body set), then fills the BW with the highest J_div; Hybrid Search (HS)
  fills half the BW by diversity and half by lowest delta-v. HS "does not search for solutions that are Pareto optimal".
  Candidate next bodies are pre-filtered with the Improved Orbital Indicator (Hennes/Izzo) before Lambert arcs are computed.
- Stated limits: the method "does not guarantee to find local optima or good Dv solutions"; pairwise distances scale badly
  with population size; a multi-objective version is deferred to future work.

## Populations and results (Q2)
Both cases: 500 tours of 10 asteroids from Earth, 5-year launch window from 1 Jan 2025, 10-day launch steps, 100-400 day legs, BW = 500.
- Main belt (16,256 asteroids, GTOC7 data): BS 2 unique tours and 10 unique asteroids, mean dv 18.87 km/s; DS 500 unique, 717
  asteroids, 42.54 km/s; HS 500 unique, 603 asteroids, 30.88 km/s and the lowest single tour (17.23 km/s vs BS 18.72).
  Diversity score mean: BS 0, DS 31.42, HS 24.93 km/s.
- Near-Earth (1,436 asteroids, GTOC4 data): BS 30 unique tours, 43 asteroids, 37.24 km/s; DS 500 unique, 376 asteroids, 50.67
  km/s; HS 339 asteroids, 46.52 km/s (HS did not beat BS on minimum dv here: 37.48 vs 33.20). Amors+Apollos visited:
  BS 43, DS 376, HS 339; no Atens or Atiras reached.

## Dynamics and legs (Q3)
Each leg is a Lambert arc between heliocentric orbits, impulsive rendezvous (match asteroid position and velocity), delta-v
limits per departure/arrival (2 km/s main belt, 3.5 km/s NEA), epochs on a discrete grid. "The current dynamical model
considers impulsive rendezvous maneuvers with asteroids where trajectories are computed with Lambert arcs". No gravity
assists in the experiments, no low-thrust, no multi-body dynamics; refinement to deep-space manoeuvres or low-thrust is left to
future work.

## Gravity assists / planetary or moon tours (Q4)
One sentence, in the discussion, as an untested claim: "the DS and HS algorithms are versatile and can be applied to various
multi-target missions, including on-orbit servicing and multiple flyby missions, also including gravity assist with Solar
System planets as in Bellome et al. (2024), as well as adaptable to any model, from multiple deep-space maneuvers to
low-thrust trajectories." Moon tours: not stated. (Observation: the reference list has Bellome et al. 2023, the gravity-assist
dynamic-programming paper, and Bellome et al. 2024, the asteroid-belt one; the citation year in this sentence matches the
latter.) The introduction also notes that dynamic programming over tree-like spaces appears in Bellome 2024 for global
optimality, as an alternative to heuristic pruning.

## Word search (case-insensitive, line breaks joined, ligatures normalised)
cycler: 0. periodic: 0. resonan: 0. three-body: 0. gravity assist / gravity-assist: 4 (related work: ant-colony search of
multi-gravity-assist transfers; the discussion sentence quoted above; reference titles Bellome 2023 and Ceriotti-Vasile 2010a).
"MGA" appears in related work only. "flyby" 5, "fly-by" 1 (related work and reference titles).

## Relevance to this project
Low. Heliocentric Lambert rendezvous chains between asteroids, no resonance, periodicity, cycler or gravity-assist
computation. Two ideas could transfer to a search over sequences of planetary or moon encounters, the paper supporting them
only as a combinatorial-search technique: (1) a set-relative diversity score used as the beam-pruning criterion, so that a
beam over encounter sequences does not collapse onto many epoch-variants of one body sequence (the paper's BS returned 2
unique tours out of 500); (2) the hybrid rule (half the beam by cost, half by diversity) as a cheap guard against a
cost-ranked search missing structurally different families, which could also be used to report coverage of a sweep. The paper
makes no claim of either for gravity-assist tours beyond the single sentence quoted.

## Corpus index line
Grabowski, Bellome & Felicetti 2026, ASR 77 7039: diversity score and beam-search variants (DS, HS) for main-belt and NEA
Lambert-arc asteroid tours; no gravity assists tested (one untested applicability sentence), no cyclers; relevance low.
