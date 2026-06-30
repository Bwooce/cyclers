# Dellnitz et al. 2005 — Transport of Mars-crossing asteroids from quasi-Hilda region (digest)

Digested 2026-06-30 (Task #499 heteroclinic/transport pass).

**Source (cite exactly):**
Dellnitz, M., Junge, O., Lo, M. W., Marsden, J. E., Padberg, K., Preis, R.,
Ross, S. D. and Thiere, B.,
"Transport of Mars-Crossing Asteroids from the Quasi-Hilda Region,"
*Physical Review Letters* **94** (2005) 231102.
DOI: 10.1103/PhysRevLett.94.231102. Received 30 September 2004; published 16 June 2005.

Full text-layer PDF read. 4 pages (PRL letter); complete text extracted.

---

## 1. One-line contribution

First application of **set-oriented (GAIO) methods** with **graph partitioning
(PARTY)** to compute *quantitative transport rates* in a celestial mechanics problem:
the fraction of quasi-Hilda asteroids (near the 3:2 Jupiter resonance) that transit to
Mars-crossing orbits, calculated without long-time single-trajectory simulations.

---

## 2. Model and method

**Model:** Planar circular restricted three-body problem (PCRTBP) with Sun and Jupiter
as primaries. Hamiltonian:
H = E = ½(p²_x + p²_y) − xp_y + yp_x − m_S/r_S − m_J/r_J − ½m_S m_J
(Eq.1; rotating frame with standard normalization: m_S = 1−µ, m_J = µ, separation = 1).
Coordinates: (x,y) position; (p_x = ẋ−y, p_y = ẏ+x) canonical momenta.

**Poincaré surface of section (SOS):**
Section defined by y = 0, ẏ < 0, x < 0 (interior realm: mean motion greater than
Jupiter's, conjunctions on opposite side of Sun from Jupiter). Coordinates on section:
(x, ẋ = dx/dt). Area- and orientation-preserving return map f: M → M, where M ⊂ ℝ².

**Mars-crossing line:** Locus of constant periapse r_p = a(1−e) = a_Mars (semimajor
axis of Mars' orbit). Points to the right of this line in (x, ẋ) space have periapsis
inside Mars's orbit → Mars-crosser.

**Energy:** E = −1.52 (just below L1 energy E_L1 = −1.5199). This choice:
- L1 neck between interior and Jupiter regions is **closed** (no dynamical connection
  between interior and Jupiter regions at this energy)
- Jupiter perturbation still significant (not zero)
- Contains the 3:2 Hilda resonance island ("U"-shaped region at left of SOS, Fig.1)

**Set-oriented method (GAIO framework):**

Step 1 — Box covering: cover the recurrent set in M with a collection of boxes B₁,...,B_b.
Choose n points per box; iterate each point forward one application of f; record which
box the image lands in. This builds a transition graph with boxes as vertices.

Transition probability from box B_i to box B_j:
p_ij = m(f⁻¹(B_i) ∩ B_j) / m(B_j)
→ sparse matrix P = (p_ij).

Step 2 — Almost-invariant set decomposition: partition the graph vertices into a small
number of parts such that each part is highly intra-coupled and loosely inter-coupled.
The key computational insight: this is equivalent to a **graph partitioning** problem
(NP-complete in general); solved with PARTY software using multi-level heuristics
(coarsen large graph → partition small graph → project back, improving at each level).
Heuristics: maximum-weighted-edge matching (coarsening) + helpful-set heuristic
(local improvement). Both provide provable quality guarantees.

The result: three almost-invariant sets found for this problem (Fig.2). One of them
contains the Hilda resonance island → this is the "quasi-Hilda region" R.

Step 3 — Transport probability: transition matrix P gives:
p_{R→Q}(n) = (1/m(R)) · e_Q^T · P^n · u_R
where e_Q ∈ ℝ^b has 1 if B_i ⊆ Q (0 otherwise), u_R ∈ ℝ^b has m(B_i) if B_i ⊆ R
(0 otherwise), and Q is the region beyond the Mars-crossing line.
This gives the probability that a particle starting in R (quasi-Hilda region) maps
into Q (Mars-crossing region) after n iterations of f.

**Tools named:**
- GAIO (Global Analysis of Invariant Objects): dynamical systems software library
  implementing the set-oriented methods [ref 10, 11]
- PARTY: graph partitioning software library [ref 15, Preis PhD thesis, Paderborn 2000]
- GADS (Graph Algorithms for Dynamical Systems): new tool that combines GAIO + PARTY
  + further graph algorithm implementations [ref 16, "to be published" at time of paper]

**Advantage over long-time simulations:**
- Classical approach: propagate many test particles for long times, extract statistics.
  Problem: chaotic dynamics + "trapping" (almost-invariant sets) → trajectories can
  be trapped exponentially long → statistics unreliable with any finite simulation set.
- GAIO approach: uses *short-time* simulations for *many globally distributed* initial
  conditions; tracks evolution of *densities* not trajectories; a priori error estimates
  available [ref 12]. Not affected by near-ergodicity failures.

---

## 3. Sourced goldens (exact numbers from the paper)

**Energy level used:**
- E = **−1.52** (just below L1: E_L1 = −1.5199)
- Caption: "The energy we consider in this Letter, E = −1.52, is just below that of
  the equilibrium point L1 (E = −1.5199)"

**Transport results (Fig.4, p.231102-3):**
- Probability for a typical quasi-Hilda particle to leave the quasi-Hilda region:
  **~6% after 200 iterates** of the Poincaré map f
- Corresponding transit time to Mars-crossing: **2000–6000 Earth years** (depending
  on location within quasi-Hilda region; the map f is crossed every Jupiter conjunction
  ≈ Jovian synodic period ≈ 398 days, so 200 iterates ≈ 218 years; the paper states
  "between 2000 and 6000 Earth years" which suggests each iterate corresponds to
  ~10–30 Jupiter conjunctions on average, i.e. multi-revolution tracking)
  NOTE: The text says "200 iterates ... corresponds to a transit time between 2000 and
  6000 Earth years" — this is the dominant timescale, not 200 × conjunction interval.

**Almost-invariant sets:** **3 sets** found for the quasi-Hilda decomposition (Fig.2).

**Computation time:** **~3.5 hours** on Intel Xeon 3.2 GHz processor (covering +
transition matrix computation dominates total cost).

**Key formula (Eq in text):**
p_{R,Q}(n) ≈ (1/m(R)) · e_Q^T · P^n · u_R
(transition probability from quasi-Hilda region R to Mars-crossing region Q after n
iterates; P is the b×b transition matrix, b = number of covering boxes)

**Hilda 3:2 resonance island:** Visible as a "sideways U"-shaped (horseshoe-shaped)
region on the left of the SOS at the plotted energy E = −1.52 (Fig.1).

---

## 4. Reuse verdicts

### #308 — Asteroid-leveraging transport (set-oriented transport rates)

**DIRECT MATCH for method.** This paper is exactly the proof-of-concept for
set-oriented transport in a celestial-mechanics context. The two contributions relevant
to #308 are:

(a) **Method availability:** GAIO + GADS software exists and has been applied to a
solar system transport problem. The quasi-Hilda-to-Mars computation is structurally
analogous to a "how quickly do unstable-resonance objects reach Earth-crossing orbits"
computation. If #308 involves estimating asteroid-leveraged assist rates or flux
into cycler-relevant orbits, this is the right computational paradigm.

(b) **The "trapping problem" warning:** Near-ergodic phase spaces (like the ones with
large KAM tori surrounding the 3:2 resonance) cannot be adequately sampled by a few
long simulations. Any transport rate claim for #308 based on Monte Carlo propagation
of a small ensemble should cite this paper as the reason for caution and GAIO as the
more reliable alternative.

**Limitations for direct reuse:**
- Model is *planar* CR3BP (authors note extension to 3D is future work)
- Only one energy level (E = −1.52) studied; authors note multi-energy robustness
  as future work
- GAIO/GADS software is academic (Paderborn group); not SPICE-compatible

### #314, #411 — Cross-system heteroclinic, 3D spatial

**NOT APPLICABLE.** This paper is set-oriented transport, not manifold-tube
trajectory design. The reference to the Lo–Ross manifold network (refs 1–3 in the
paper) is background motivation only. No trajectory construction, no manifold tubes,
no 3D extension.

### #291, #306 — 3D families, spatial manifolds

**NOT APPLICABLE.** Planar method; no manifold computation.

---

## 5. Key algorithmic detail for future implementation

The GAIO box-covering method in three practical steps:
1. Choose box size and number of sample points per box. Grid the recurrent set M.
2. For each box: sample n_pts points uniformly, apply f once, record destination box.
   This gives the transition graph (b vertices, weighted directed edges).
3. Partition the graph using PARTY (or any multilevel graph partitioner): find a
   k-way partition (k = number of almost-invariant sets desired) minimizing the
   edge-cut weight between parts.

The transition matrix P is then used directly for transport probability calculations
via matrix power P^n (or iterative matrix-vector products for large b).

**Key parameter:** box size δ and sample count n_pts control the accuracy/cost tradeoff.
A priori error estimates (ref [12] = Dellnitz–Junge–Koon–Lekien–Lo–Marsden–Padberg–
Preis–Ross–Thiere, *Int. J. Bifurcation Chaos* **15** (2005) 699) bound the
approximation error.

---

## 6. Status

**Digested** (full paper read, all sourced numbers extracted). CORPUS_INDEX update
deferred to the calling task's commit step.
