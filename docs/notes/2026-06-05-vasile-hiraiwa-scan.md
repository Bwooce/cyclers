# Relevance scan: Vasile & Campagnola (2009) + Hiraiwa et al. (2026)

Date: 2026-06-05. Scan-depth only (not deep mining). Per-paper verdicts below.
Citations by author/title/arXiv ID only; PDF paths deliberately omitted (private mirror).

---

## Paper 1 — Vasile & Campagnola, "Design of Low-Thrust Multi-Gravity Assist Trajectories to Europa" (JBIS; arXiv:1105.1823)

PDF has broken font encoding (pdftotext garbage); read via 150/220-dpi page rasters.
Pages read: title/abstract (p.1-2), results figs + Tables 3 & 4 + Conclusions (p.22-25). Did NOT read the
method body (p.3-21) in detail — only abstract-level understanding of the algorithm.

### Cyclers?
No periodic two-body-encounter cycler orbit. This is a one-shot interplanetary delivery trajectory
(Earth→…→Jupiter→Europa). However it contains a **planet-centric multi-moon gravity-assist tour** of the
Jovian system that is directly relevant to deferred interest #76. Abstract (p.1): "an additional manoeuvre
would be use solar electric propulsion or a combination of impulsive and chemical propulsion." Figure
captions explicitly name a **"synchronous tour of Ganymede"** (Fig.17) and a **"synchronous tour of Europa"**
(Fig.18) — i.e. resonant moon-flyby sequences, the moon-tour analogue of a cycler, not a closed cycler.

### Numeric results (candidate golden anchors for low-thrust MGA, task #37)
Full optimised solution is an MGA-LT chain. Body sequence (Table 3, "Time of Flights and encounter dates
for the entire optimised solution"; dates in MJD, TOF in days, cumulative mission time in years):

| Ph | Body          | Dep date (MJD) | TOF (d) | Cum. (yr) |
|----|---------------|----------------|---------|-----------|
| 1  | Earth         | 3718.9         | 393.7   | 1.08 |
| 2  | Venus         | 4112.7         | 166.5   | 1.53 |
| 3  | Mars          | 4279.1         | 375.3   | 2.56 |
| 4  | Earth         | 4634.4         | 1241.6  | 5.96 |
| 5  | SOI of Jupiter| 5896.1         | 110.6   | 6.27 |
| 6  | Ganymede1     | 6006.7         | 492.7   | 7.62 |
| 7  | Ganymede2     | 6499.3         | 93.0    | 7.87 |
| 8  | Ganymede3     | 6592.5         | 35.8    | 7.97 |
| 9  | Ganymede4     | 6628.1         | 21.3    | 8.03 |
| 10 | Ganymede5     | 6649.6         | 1.2     | 8.03 |
| 11 | Europa1       | 6650.8         | 10.7    | 8.06 |
| 12 | Europa2       | 6661.4         | 17.8    | 8.11 |
| 13 | Europa3       | 6679.2         | 7.1     | 8.13 |
| 14 | Europa4       | 6686.3         | 28.4    | 8.21 |
| 15 | Europa5       | 6714.7         | 39.6    | 8.32 |
| 16 | Ganymede6     | 6754.4         | 2.0     | 8.32 |

(Reading at raster resolution; treat MJD/TOF digits as approximate — verify against the published JBIS table
before using as a hard golden value.)

Table 4 ("Summary of gravity assist characteristics for the entire optimised solution") gives per-flyby
departure/arrival relative & absolute velocities, bend angle B [deg], pericentre altitude, and Δv_GA [km/s]
for the same 16 phases. Representative rows: Earth→Venus dep rel-V 1.78, arr rel-V 0.71 km/s, B=36.89°,
peri-alt 753 km, Δv_GA 6.88; Venus→Mars 0.71/0.71, B=20.08, alt 17478, Δv 0.60; the deep-space Ganymede tour
rows (phases 6-16) carry Δv_GA in the 0.5-1.1 km/s range with pericentre altitudes ~200-580 km.

### Search algorithm (one paragraph)
Two-layer hybrid global optimisation. An outer **genetic algorithm (GA)** searches the discrete/combinatorial
space of flyby phasing and resonant gravity-assist sequences; the inner layer transcribes each leg's
low-thrust optimal-control problem with **Direct Finite Elements in Time (DFET)** and solves it
deterministically. Abstract (p.2): the GA "generates families of resonant gravity assist trajectories" and the
DFET transcription "presents a number of first guess solutions and a fully optimised transfer for Europa."
This GA-branches-then-deterministic-refine pattern is the relevant idea for a trajectory-search pipeline.

### Verdict — Paper 1
- Catalogue rows: **not a cycler**, so no direct cycler-catalogue row. But the Ganymede/Europa **synchronous
  (resonant) moon tour** is a strong reference for deferred interest #76 (planet-centric moon-tour cyclers) —
  worth a future "near-cycler / moon-tour" annotation, not a periodic-cycler row.
- Golden anchors: **possible** for low-thrust MGA (#37) from Tables 3 & 4, but ONLY if the numbers are
  re-read from a clean copy of the JBIS table (raster digits are not trustworthy as goldens). Flag as
  "needs clean-source transcription" before any golden test uses them.
- Algorithm ideas: yes — GA over resonant-flyby sequences + DFET inner solve.

---

## Paper 2 — Hiraiwa, Bando, Sato, Hokamoto, "Design of low-energy transfers in cislunar space using sequences of lobe dynamics" (arXiv:2602.17444, 2026; submitted to Acta Astronautica)

Normal text PDF. Pages read: highlights/abstract/intro (p.1-2), Section 6 results incl. Tables 2 & 3 and
Figs 28-32 (p.33-40), conclusions + references (p.40-47). Did NOT read the lobe-dynamics theory (p.3-32) closely.

### Cyclers?
No. This is a single LEO→LLO Earth–Moon low-energy transfer in the CR3BP/BCR4BP — not a periodic
multi-encounter cycler. Relevant instead to the existing Earth–Moon CR3BP machinery (deferred CR3BP rows).
Abstract (p.1): "The resulting optimal trajectory in the Earth–Moon CR3BP is then converted into an optimal
transfer in the bicircular restricted four-body problem (BCR4BP) via multiple shooting."

### Numeric results (candidate golden anchors — Earth–Moon transfers / CR3BP–BCR4BP)
Table 2 (p.34, "Results of the optimization in the BCR4BP"; LEO 167 km → LLO 100 km; θ_s* = Sun initial phase):

| θ_s*  | Total ΔV [m/s] | Transfer time [day] |
|-------|----------------|---------------------|
| 0     | 3832.6088      | 193.2512 |
| π/6   | – (no conv.)   | – |
| π/3   | –              | – |
| π/2   | 5201.0019      | 191.7576 |
| 2π/3  | 3841.6832      | 193.7585 |
| 5π/6  | –              | – |
| π     | –              | – |
| 7π/6  | –              | – |
| 4π/3  | –              | – |
| 3π/2  | 4934.6452      | 192.5052 |
| 5π/3  | 3848.4479      | 193.9369 |
| 11π/6 | 3837.8299      | 193.3733 |

Best solution θ_s*=0: total ΔV 3832.6088 m/s, 193.25 d. CR3BP initial-guess transfer time 191.8994 d.
For θ_s*=0 the two dominant burns are |ΔV_E| = 3120.8662 m/s (LEO departure) and |ΔV_M| = 638.5176 m/s
(LLO arrival); remaining trim burns each |ΔV| < 21 m/s (p.35).

Table 3 (p.36, "List of known Earth–Moon interior transfers, LEO 167 km → LLO 100 km" — literature
comparison, sourced from Topputo 2013 [22]). These trace to published prior work, so they are usable as
external golden/comparison anchors:

| Reference                     | Total ΔV [m/s] | Time [day] |
|-------------------------------|----------------|-----------|
| Hohmann transfer [22]         | 3954           | 5 |
| Sweetser [100] (theoretical min)| 3726         | – |
| Pernicka et al. [101]         | 3824           | 292 |
| Yagasaki [102]                | 3925 / 3947 / 3951 | 31 / 14 / 4 |
| Yagasaki [103]                | 3941 / 3949    | 14 / 5 |
| Topputo et al. [104]          | 3895 / 3900    | 256 / 194 |
| Mengali & Quarta [105]        | 3861 / 3920 / 3950 / 4005 | 85 / 68 / 14 / 3 |
| Mingotti et al. [106]         | 3896 / 3917 / 3936 | 31 / 30 / 14 |
| Da Silva Fernandes & Marinho [107] | 3850 / 3902 / 3943 / 3950 | 58 / 32 / 14 / 5 |
| Topputo [22]                  | 3893 / 3937 / 3945 | 31 / 14 / 5 |

CR3BP energy markers along the optimal trajectory: L1 Lyapunov orbit Jacobi constant C_J = 3.16; perturbed
C_J along the BCR4BP trajectory stays in 3.1551 < C_J < 3.1666 (p.37-39). Resonance ratios crossed by the
semi-major-axis history: 7:2, 3:1, 5:2 (Fig.30).

### Search algorithm (one paragraph)
**Lobe-dynamics graph search.** Rather than detecting full lobe transport structures numerically, the method
introduces "effective lobes" as robust intermediate transfer points between start and goal orbits. A
weighted, directed graph is built whose nodes are effective lobes and whose edges are transfer arcs
(determined by a targeting strategy); a combinatorial optimisation over this graph selects a sequence of
lobes that connects departure and arrival orbits at low energy (Conclusions, p.40). The chaotic CR3BP
trajectory threaded through the chosen lobes is then refined into a BCR4BP solution by **multiple shooting**
(fmincon / SQP, N=45 segments, 12 Sun-phase initial guesses θ_s*=kπ/6; p.33). This lobe-sequencing-as-graph
idea is the transferable algorithm for a trajectory-search pipeline.

### Verdict — Paper 2
- Catalogue rows: **no cycler rows** (single low-energy transfer, not periodic/multi-encounter).
- Golden anchors: **yes, good ones** for Earth–Moon transfer / CR3BP–BCR4BP machinery. Table 3 values are
  literature-sourced (independent of this paper's own code) → suitable as golden EXPECTED values per the
  "sourced-only" rule. Table 2 (this paper's own optimiser output) is usable as a regression/reference anchor
  but is self-computed, so NOT a golden EXPECTED value. C_J=3.16 (L1 Lyapunov) is a clean CR3BP anchor.
- Algorithm ideas: yes — effective-lobe directed-graph sequencing + multiple-shooting CR3BP→BCR4BP refine.
