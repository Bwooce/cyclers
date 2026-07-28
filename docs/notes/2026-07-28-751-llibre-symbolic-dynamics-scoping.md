# #751 — Llibre-Martínez-Simó 1985 Theorem C vs. the binary-cycler genome: scoping verdict (2026-07-28)

Research/scoping only (no code, no catalogue writes, no runs), mirroring `#714`'s
format: an honest tractability/worth verdict first, mechanics second. The question
(user-flagged after the `#749` digest): does Theorem C's Bernoulli-shift
symbolic-dynamics existence proof of quasirandom Sun<->Jupiter shuttling orbits
represent a genuinely different construction method from this project's existing
"binary-cycler genome" tooling, or the same ground under different formalism — and
is it worth reopening capability-building
(`[[project_capability_frontier_complete]]`) for?

Inputs actually read (not summarized from memory): the Llibre-Martínez-Simó 1985
paper itself, §1 (Theorem C's simplified statement + context, pp. 108-110) and the
FULL §6 (pp. 147-155: the two-big-squares/Devaney singular-cross-section setup,
Proposition 6.1, Theorem 6.1 + proof, Lemma 6.1 + proof, Theorem C full statement +
proof, Theorem 6.2, and the closing periodic-family remark), plus §5's Table I
(`mu_k` values) and the (46,46)-homoclinic worked point (`mu=0.00025,
Delta C=0.0000783`), from
`/Users/bruce/dev/cyclers_pdf/papers/llibre-martinez-simo-1985-transversality-invariant-manifolds-lyapunov-l2-jde-58-104-doi-10.1016-0022-0396(85)90024-5.pdf`
(via its `.txt` sidecar); the `#749` digest in full; the genome package itself
(`genome/__init__.py`'s genome definition, `search/binary_star_search.py`'s
`Topology`/`winding_topology` in full, `genome/asymmetric_branch.py` header,
`genome/heteroclinic_cycle.py` header + `symbol_sequence` machinery,
`genome/transit_manifold.py` header, `genome/composed_moon_map.py`'s
`run_itinerary`); the `#746` BMO/Casoliva-He1 cross-check note in full; and the
`#714` scoping-pass note as the format precedent.

## 1. What Theorem C actually gives (read from the proof, not the abstract)

**Statement (paper's own final form, §6 p. 154).** Fix `(mu, C)` such that the two
branches `W^u_+, W^u_-` of the L2 Lyapunov orbit's unstable manifold meet the
corresponding stable branches transversally (L2 in the paper's Szebehely
convention = the collinear point BETWEEN Sun and Jupiter; `L` = the neck region,
`S`/`J` = the Sun-side/Jupiter-side realms). Then for ANY doubly-infinite integer
sequence `a = (..., m_-1, m_0, m_1, ...)` with `|m_j| > m = m(mu, C)` for all `j`,
there exist initial conditions — **unique in a neighborhood of the chosen
homoclinic orbits** — whose orbit realizes the itinerary: on pass `j` it enters
`L`, makes exactly `|m_j|` revolutions around the Lyapunov orbit, and exits to the
same side or crosses to the opposite side according to `sign(m_j)`. The
compactified sequence types beta/gamma/delta additionally give orbits forward-,
backward-, or bi-asymptotic to the Lyapunov orbit itself.

**Answers to the dispatch's specific sub-questions:**

- **Constructive or pure existence?** Pure topological existence + local
  uniqueness. The proof instantiates Moser's strip machinery (Theorem 6.1): the
  orbit is `h(a)` = the single intersection point of one vertical and one
  horizontal curve, each obtained as the limit of NESTED strips
  (`V = F^{-1}_{|a_1|}(A_{i_1} ∩ F^{-1}_{|a_2|}(...))`) whose diameters contract
  exponentially (condition (b): `d(V') <= s·d(V)`, `s < 1/2`, verified via the
  sector-field condition (c) which holds because the neck flow is EXACTLY linear
  in Moser-Lyapunov coordinates). There is **no corrector, no continuation
  recipe, no error bound usable as an algorithm** — the strips are defined
  through infinitely many iterates of the (chaotic) return map, which is
  precisely what a numerical method cannot iterate naively. Any code
  implementation must invent its own numerical construction (the standard route:
  finite-block itinerary -> segment chain along numerically-globalized
  manifolds -> multiple-shooting shadowing correction). The local-uniqueness
  clause is numerically GOOD news (well-posed shadowing target), but the theorem
  contributes the guarantee, not the method.
- **Parameter regime.** The rigorous chain is asymptotic: Theorem A gives the
  countable homoclinic sequence `mu_k = 1/(N(inf)^3 k^3)(1+o(1))` as `mu -> 0`;
  Theorem B gives transversal homoclinics for
  `Delta C > L·mu_k^{4/3}(mu - mu_k)^2` NEAR each `mu_k`; Theorem C needs only
  "a transversal homoclinic pair exists at `(mu, C)`". Table I (read directly):
  `mu_2 = 4.2539e-3`, `mu_3 = 6.7525e-4`, `mu_4 = 2.1929e-4`, decreasing ~`k^-3`
  to `mu_23 = 6.64e-7`. **The physical Sun-Jupiter mass ratio (9.5388e-4, this
  project's `cr3bp_system("Sun","Jupiter").mu`) sits BETWEEN `mu_3` and `mu_2`,
  at none of the `mu_k`** — at the real mass ratio, Theorem C's applicability
  rests on the paper's NUMERICS, not its theorems: Table II computes the
  first-tangency boundary `Delta C(mu)` over `mu ∈ [0.00023, 0.0042]` (which
  brackets the real value), and §5's worked point (`mu = 0.00025,
  Delta C = 0.0000783`, (46,46)-homoclinic points in the J region, "For those
  values of mu, Delta C Theorem C will be applicable") is the paper's own
  positive control. The regime is `C` just below `C2` (neck barely open) — the
  slow, high-energy-precision bottleneck regime.
- **A limitation the digest did not surface, load-bearing for cycler design:**
  the theorem's alphabet controls ONLY (side-switching sign, neck revolution
  count). The outer-region behavior is FIXED: "the number of integer revolutions
  of the massless body around the Sun or Jupiter... is a constant depending on
  the region and on the selected homoclinic orbits." **You cannot prescribe the
  S-region or J-region resonance sequence through Theorem C's symbols** — the
  outer resonance is inherited once-and-for-all from the chosen homoclinic pair.
  A KLMR-style resonance-transition itinerary (which resonance on which pass) is
  exactly what this alphabet does NOT index.
- **The cycler-relevant corollary** is the paper's closing remark: periodic
  symbol sequences give "a countable set of families of simple periodic orbits
  near the homoclinic one" — i.e., homoclinic-shadowing periodic orbits. That
  class ALREADY has a constructive, numerical, published instantiation in this
  project's corpus: **Casoliva et al. 2008/2010's Class 2 (He1-4/Hm1-2)
  Earth-Moon cyclers** (`#725`), whose theoretical-existence underpinning the
  `#749` digest already identified as exactly this paper's machinery.

## 2. Same ground as the binary-cycler genome? NO — different phase space, different role for the integers

Judged against the actual code, not the package name:

- **What the genome actually is** (`genome/__init__.py`, read directly): "compact,
  parametric descriptions of orbit families where the family member is selected
  by a small set of integer / continuous labels." The "binary-cycler genome"
  concretely = `search/binary_star_search.py`'s `Topology` dataclass: winding
  numbers `(k1, k2)` about each primary (+ `k_z` plane-crossings for 3D),
  computed by integrating a FOUND periodic orbit and unwrapping `arctan2`. The
  search itself is over ICs (figure-read seeds, correctors,
  pseudo-arclength/`mu` continuation); the integer labels are an **a-posteriori
  classification invariant** used as an independent cross-check gate
  (`[[feedback_orbit_closure_discipline]]`), not a search index. The orbit class
  is Roberts-Tsoukkas-Ross STABLE prograde multi-orbiter cyclers — regular,
  non-hyperbolic phase space, nowhere near the L2 neck's homoclinic tangle.
- **Theorem C's integers are the opposite object**: an a-priori COMPLETE index —
  every admissible bi-infinite sequence is realized by a (locally unique) orbit,
  before any orbit is computed. The orbits are chaotic transit orbits (or, for
  periodic sequences, unstable homoclinic-shadowing periodic orbits) living
  precisely in the hyperbolic tangle at `C` just below `C2`.
- **Closest existing code is NOT the binary genome but the #314/#547 tube-dynamics
  modules**, and neither covers this ground either:
  `genome/heteroclinic_cycle.py` certifies finite closed HETEROCLINIC cycles
  between DISTINCT L1/L2 Lyapunov orbits (its `symbol_sequence` is a
  certification OUTPUT for a found cycle, not a prescribed-itinerary
  constructor; its own docstring's "symbolic dynamics" reference is the
  unrelated Wilczak-Zgliczyński L1<->L2 result, as `#749` already noted);
  `genome/transit_manifold.py` classifies transit vs non-transit branches
  (no homoclinic return at all); `genome/composed_moon_map.py`'s
  `run_itinerary` propagates Keplerian-map encounter itineraries — a
  patched-map, not manifold-tangle, construct. And `#746` independently
  confirmed the project has **no Lyapunov-orbit homoclinic-connection
  finder/continuator of any kind** (the `ccr4bp_heteroclinic_search.py`
  least-squares refiner is architecturally different and heteroclinic).

**Verdict (a): genuinely different ground.** Bi-infinite-itinerary-indexed
orbits near the L2-Lyapunov homoclinic tangle are not searched, encoded, or
certified anywhere in this codebase. The overlap with the binary-cycler genome
is superficial (both use integer labels); the roles, orbit classes, and
phase-space regions are disjoint.

## 3. What a build would take, and what the theorem contributes to it

Because Theorem C is non-constructive, a "symbolic-dynamics-indexed orbit
constructor" is really the standard numerical itinerary-shadowing pipeline —
which is **published prior art, constructively executed by others**: KLMR 2000
(Sun-Jupiter resonance transitions with prescribed `(J,X;J,S,...)` itineraries,
in corpus) and Casoliva 2008/2010 (Earth-Moon periodic homoclinic-shadowing
cyclers, in corpus, with golden numeric Tables 4-6 — flagged by `#725`/`#746` as
"the first sourced numeric target" for any future homoclinic build). Concretely:

1. **Homoclinic-connection finder** for a single Lyapunov orbit: globalize
   `W^u`/`W^s` at fixed `C` (Floquet seeding exists in `heteroclinic_cycle.py`),
   intersect on the `y=0` section, refine the transversal homoclinic point.
   Missing today; moderate. **~3-5 days** including the Casoliva Table 4 golden
   reproduction (He1, `h = -1.45016232260699`, 113.6-day connection, Earth-Moon
   L1) and Llibre's own §5 controls (`mu_k` Table I; the (46,46) point).
2. **Itinerary -> orbit shadowing constructor** for finite/periodic symbol
   blocks: chain homoclinic-loop segments with `m_j` neck revolutions inserted
   per pass, close with multiple shooting. `cr3bp_multiple_shooting.py` /
   `genome/multi_shooting.py` exist but are shaped for periodic-orbit closure
   (`#746`'s node/segment-time contract caveat) — adaptation, not reuse.
   **~1-1.5 weeks**, dominated by conditioning in the near-`C2` regime (segment
   times through the neck grow like the theorem's own `m` lower bound).
3. **Certification + ratchets**: extend `winding_topology`-style itinerary
   verification per pass, periodicity residual for cyclic sequences, ghost-guard
   style independent-integrator check. **~2-3 days.**

**Verdict (b): ~2-3 weeks total** for the full constructor; a minimum viable
slice (item 1 alone: homoclinic finder + Casoliva He1 golden reproduction) is
**~1 week**. Same order as the `#714` N=5 estimate, but with a materially worse
payoff profile (below).

What the theorem DOES contribute if a build ever happens (the "even if same
ground" question, answered even though the ground is different): (i) a rigorous
realizability guarantee — every sequence with `|m_j| > m(mu,C)` exists and is
locally unique, so a shadowing failure above that floor is a solver bug, not
physics (a qualitatively new kind of gate for this project: a
literature-sourced EXISTENCE lower bound rather than a numeric golden value);
(ii) the `mu_k`/`N(inf)`/`M(inf)`/Table I-II numbers as digit-grade positive
controls for step 1; (iii) the `m(mu,C)` floor and fixed-outer-resonance
constant as free consistency checks on any constructed member.

## 4. GO/NO-GO

**NO-GO on registering a Theorem-C-driven symbolic-dynamics-constructor
capability build as a standalone task.** Reasoning, in order of weight:

1. **The discovery payoff is reproduction-shaped, not discovery-shaped.** The
   constructive version of this theorem's content is already published: KLMR
   2000 (aperiodic/transit itineraries, Sun-Jupiter) and Casoliva 2008/2010
   (periodic homoclinic-shadowing cyclers, Earth-Moon). The first outputs of a
   2-3-week build would be reproductions; "novel" outputs would be
   itinerary-level enumeration deeper into the same published families (longer
   symbol blocks, other planet-moon pairs) — weak novelty of exactly the kind
   the `[[feedback_literature_novelty_check_baseline]]` gate exists to
   deflate. Contrast `#714`'s GO: there, the literature itself said no
   periodic orbits were known in the target model (genuinely open ground), and
   it produced a banked novel discovery.
2. **The alphabet doesn't index what cycler design needs.** Theorem C's symbols
   control neck revolutions + side-switching only; the outer-region resonance —
   the quantity this project's actual cycler encodings (S/L resonance
   intervals, encounter sequences, `(k1,k2)` windings) are built around — is
   frozen by the choice of homoclinic pair. As a cycler-search index it is
   strictly less expressive than what already exists; its strength
   (completeness guarantee) is a validation asset, not a search axis.
3. **The regime is operationally marginal**: `C` just below `C2`, revolutions
   floor `m(mu,C)`, TOFs long and chaotic — quasi_cycler-class objects at
   best, in the hardest-to-verify corner of phase space.

**One narrow, deferred alternative — do NOT open it independently of `#752`:**
the genuinely missing primitive both this scoping AND `#746` AND `#742`
independently identified is a **Lyapunov-orbit homoclinic-connection
finder/continuator** (step 1 above, ~1 week, with Casoliva Tables 4-6 as sourced
golden targets and Llibre Table I/II as regime controls). That primitive is
shared substrate with `#752`'s Anderson/Lo resonant-orbit-manifold/homoclinic
tour-design scoping (running concurrently). If `#752` comes back GO, fold the
homoclinic-connection finder into that build and register Llibre's constants as
its positive-control battery; if `#752` also comes back NO-GO, this paper's
value is fully banked as-is (digest `#749` + this note + the existence-guarantee
framing) at zero build cost. No new task number self-assigned here, per scope.
