# Digest - Zhang, Guo, Wu, Baoyin, Li & Topputo (2025), "Global Optimality in Multi-Flyby Asteroid Trajectory Optimization: Theory and Application Techniques"

**Digested:** 2026-10-03 (text-layer PDF; all 39 pages read). **Citation:** Zhong Zhang (Tsinghua; Politecnico di Milano),
Xiang Guo (NUDT), Di Wu (Beihang), Hexi Baoyin and Junfeng Li (Tsinghua; Li is the corresponding author), Francesco Topputo
(Politecnico di Milano), "Global Optimality in Multi-Flyby Asteroid Trajectory Optimization: Theory and Application
Techniques". Author manuscript (revision R1 per the file name; PDF created 17 October 2025, LaTeX/MiKTeX, AIAA-style layout
with "Member AIAA" footnotes). The manuscript prints no journal name, volume, date, DOI or arXiv number: all "not stated".
The publication venue is therefore not established from the file. Funding: National Natural Science Foundation of China
(grants 12372047, 62227901). Code and data: GitHub repository `zhong-zh15/Multi_Flyby_Dynamic_Programming` (stated in the
manuscript). Filed in the private paper corpus as
`zhang-guo-wu-baoyin-li-topputo-2025-global-optimality-multi-flyby-asteroid-trajectory-optimization-r1-manuscript.pdf`
(md5 `f2a59671e5a3eb304ae9a75566831af1`, 39 pages).

## What it is
A fuel-optimisation method for the second step of asteroid-tour design: given a fixed sequence of asteroid flybys (or
rendezvous), choose the flyby epochs, flyby velocities and masses to minimise propellant. Sequence selection (step one) is
explicitly out of scope: "This paper considers the second step, with a focus on designing fuel-optimal trajectories for a
given flyby sequence."

## Method (Q1)
- The original optimal-control problem P0 (flybys as interior-point equality constraints) is recast as an N-stage decision
  problem P1. State at flyby k is s_k = [t_k; v_k; m_k] (epoch, spacecraft velocity, mass); the stage cost g_k(s_k, d_k) is the
  globally optimal transfer cost between consecutive flybys. This rests on a stated Markov property: "only the velocity and
  mass of spacecraft at the flyby epoch influence subsequent trajectories".
- P2 discretises the states (time, velocity, mass grids) and replaces g_k with an approximation g~_k = g_k + eps_k (analytic
  tools or neural-network predictors). Dynamic programming (Bellman) solves P2. A stochastic variant P3 is used only for the
  proof.
- Specialisations: a bi-impulse algorithm that drops the velocity from the state (cost scales as Nt^3 in the time-step
  count), and an adaptive step refinement within a narrow tube around the previous result. The paper states the adaptive step
  "does not theoretically guarantee global optimality".
- Error bound (Eq. 42): J1[P2](s1) <= J1[P0](s1) + N * eps_max, with eps_max the largest single-stage error. Central claim,
  abstract: "the method provides a quantifiable bound on global optima errors introduced by discretization and approximation
  assumptions". eps_max is estimated empirically by random sampling (not derived analytically); the paper calls the bound
  conservative.

## Test problems and results (Q2)
Only the flyby sequence is taken from prior best solutions; epochs and velocities are re-optimised.
- GTOC4 impulsive variant (48 flybys + 1 rendezvous, bi-impulse): 25,210.4 m/s, final mass 636.7 kg, "fuel savings of 20.2 kg
  over the winning result" (Moscow State University); fixed-step-only runs stay far off even at 1-day steps (27,002.8 m/s).
- GTOC11 (ten motherships, sequences of 35-42 bodies, flyby relative speed capped at 2 km/s): improved over the winner's
  bi-impulse result on all ten, by 11.2 to 646.8 m/s per mothership, about 3.3-4.5 s CPU each. Mission 1 error bound about
  300 m/s, so "requires at least 13,640 m/s".
- GTOC4 low-thrust (University of Jena sequence, 49 flybys + 1 rendezvous): 35.9 h on 2 GPUs and a 32-core CPU, neural-network
  estimator; 19.9 kg fuel remaining, "the currently known optimal solution"; estimated error bound about 4000 m/s (conservative).

## Dynamics model (Q3)
Heliocentric two-body only (Sun gravity, Eq. 3), impulsive or low-thrust: "for the test cases presented in Sec. V, continuous
and impulsive thrust is considered in a two-body gravitational field." No planetary or lunar gravity, no third-body terms.

## Gravity assists (Q4)
Not addressed. A "flyby" here is a position match with the target body: "the flyby conditions, which require that the
spacecraft's position coincides with the flyby target at each flyby epoch without imposing restrictions on the velocity".
The target does not act on the spacecraft: velocity changes only via thrust, and there is no flyby-geometry or turn-angle
model. Phrase to remember: the general dynamics are "given in the general case without specifying the spacecraft's dynamical
environment", but nothing in the formulation or results supports a velocity-changing gravity-assist flyby. Rendezvous is
treated as a special case (adds velocity matching).

## Connection to Bellome et al. 2023 (Q5)
Cited. Introduction: "Bellome et al. [31, 32] proposed a dynamic programming method for bi-impulse solutions in multi-flyby
problems. To enable fast computation, they selected the best N solutions to enhance practical applicability, but at the expense
of theoretical global optimality." [31] is Bellome, Sanchez, Felicetti and Kemble, "Multiobjective Design of Gravity-Assist
Trajectories via Graph Transcription and Dynamic Programming", J. Spacecraft Rockets 60(5), 2023; [32] is the 2024 Acta
Astronautica "Modified Dynamic Programming for Asteroids Belt Exploration". The present paper positions itself as keeping
all states (no best-N truncation) and adding an error bound. It also extends the authors' own multi-rendezvous work [33]
(Zhang et al., JGCD 2024). The paper does not itself discuss gravity assists when citing [31].

## Word search (case-insensitive, line breaks joined, ligatures normalised)
cycler: 1 (only in reference [3], title "Asteroid Flyby Cycler Trajectory Design Using Deep Neural Networks", Ozaki et al.,
cited for neural-network approximation and DESTINY+). periodic: 0. resonan: 0. gravity assist / gravity-assist: 2, both in
reference titles ([31] Bellome 2023, and [43] "Lunar Gravity Assist in the Circular Restricted Three-Body Problem"). three-body:
2, reference titles only ([20] hybrid DDP in the CR3BP, and [43]). "flyby": 119 (the paper's own sense: asteroid position-match).
No body-text use of any of the five terms.

## Relevance to this project
Low. It is two-body, massless-asteroid flyby optimisation for a fixed given sequence; it does not search sequences, has no
gravity assists, no resonance, no periodicity and no cycler concept, and its result depends on discretised epoch grids
and (in Case 3) a trained neural-network cost estimator. Two ideas could transfer to a search over sequences of planetary
or moon encounters, but only if flyby velocity change is added, which the paper does not do: (1) the Markov-stage
formulation, where the state at each encounter (epoch, velocity or v-infinity, mass) is all the future depends on, gives an
exact stage-wise dynamic programme instead of one global NLP; (2) an explicit, checkable discretisation-error bound
(N times single-stage error), with eps_max estimated by sampling, would give a defensible "no better solution exists
within X" statement for any sweep-based negative result. Both are the authors' stated contributions; their applicability
to gravity-assist stages is our inference, not a claim in the paper.

## Corpus index line
Zhang et al. 2025 (manuscript, venue not stated): dynamic-programming global optimisation with an N*eps_max error bound for
fixed-sequence multi-asteroid flyby trajectories (two-body, impulsive and low-thrust, GTOC4/GTOC11); no gravity assists or
cyclers; cites Bellome 2023; relevance low.
