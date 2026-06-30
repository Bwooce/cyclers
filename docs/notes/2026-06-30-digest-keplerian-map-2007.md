# Ross & Scheeres 2007 — Multiple Gravity Assists, Capture, and Escape in the RTBP (SIADS 6-3)

Digested 2026-06-30. **Keystone paper for #500** (controlled Keplerian-map genome).
Also informs #318 (joint moon-tour) and #465 (multi-rev leveraging).

**Source (cite exactly, no file path):**
Ross, S.D. and Scheeres, D.J., "Multiple Gravity Assists, Capture, and Escape in the
Restricted Three-Body Problem," SIAM Journal on Applied Dynamical Systems (SIADS),
Vol. 6, No. 3, pp. 576–596, 2007. DOI: 10.1137/06065195X

> Text layer present throughout. Equations rendered partially as ASCII (fractions,
> subscripts degraded) but all numbered equations are recoverable from context. Full
> paper (20 pp. + references) read in its entirety. All figures are ASCII-only; numeric
> values are embedded in captions or text.

---

## 1. Contribution in three lines

Derives a **family of 2D symplectic twist maps** (the "Keplerian map") to approximate a
particle's motion in the Planar CR3BP at each periapsis passage, via Picard's method of
successive approximations. The map state is (ω, K): the **azimuthal angle of periapsis
in the rotating frame** (= particle–perturber separation angle at closest approach) and
the **Keplerian energy**. Maps the connected chaotic zone of the phase space, identifies
lanes of fast semimajor-axis migration, and analytically connects multiple resonant
gravity assists to ballistic capture and escape via Conley-McGehee tubes.

---

## 2. The Keplerian map: full definition

### 2.1 Hamiltonian setup (§2)

PCR3BP in inertial frame centered on m₁-m₂ barycenter, µ = m₂/(m₁+m₂) ≪ 1.
Hamiltonian (rotating frame, eq. 2.5):

    H_rot(a, e, ω̄, ν) = K(a) + µR(a, e, ω̄, ν) - G(a, e)

where K = -1/(2a) is the Keplerian energy, G = √(a(1-e²)) is the angular momentum,
ω̄ = ω - t is the argument of periapsis measured in the rotating frame.

**Jacobi constant**: C_J = -2H_rot = constant along trajectories. Tisserand parameter
approximation (valid when O(µ) terms are small):

    C_J ≈ 1/a + 2√(a(1-e²))

### 2.2 Energy kick function (§3)

The **energy kick per periapsis passage** is computed by integrating the Lagrange
planetary equations over one orbit via first-order Picard iteration:

    ΔK ≈ ΔG = -µ ∫_{-T/2}^{T/2} (∂R/∂ω) dt

where T = 2πa^{3/2} and the integral is evaluated from apoapsis to apoapsis with the
perturbation concentrated at periapsis (t* = 0).

**Energy kick function** (eq. 3.11, ignoring O(µ) correction ΔR):

    f(ω, K) = ΔK_{C_J}(ω, K) / µ ≈ ΔG / µ

Key property: **f is odd in ω**, i.e. f(-ω) = -f(ω) (antisymmetric about the m₁-m₂
line). Computed numerically on a grid of ω ∈ [-π, π] for fixed K̄; stored in a
look-up table with interpolation.

Physical meaning:
- ω > 0 (periapsis passage slightly AHEAD of moon in its orbit) → f < 0 → energy
  decreases → semimajor axis decreases
- ω < 0 (periapsis passage slightly BEHIND moon) → f > 0 → energy increases →
  semimajor axis increases
- Maximum kicks occur at ±ω_max; orbit with smaller a yields larger f_max

### 2.3 The map equations (§4, eq. 4.2) — THE GENOME

Making f independent of K (fix reference K̄) gives the area-preserving twist map
**F: (ω_n, K_n) → (ω_{n+1}, K_{n+1})**:

    ω_{n+1} = ω_n - 2π(-2K_{n+1})^{-3/2}  (mod 2π)
    K_{n+1} = K_n + µ f(ω_n)

Twist condition (right twist):

    ∂ω_{n+1}/∂K_n |_{ω_n} = 3π(-2K_{n+1})^{-5/2} > 0  for K_n < -µf(ω_n)

The map takes one periapsis passage to the next in terms of (ω, K). For given (C_J, µ),
the eccentricity e at each step follows from the Tisserand relation.

**Unperturbed map** (µ = 0, eq. 6.2): all orbits on Rotational Invariant Circles (RICs);
each point maps to (ω - 2π(-2K)^{-3/2} mod 2π, K). The KAM theorem guarantees that most
RICs persist for small µ since f integrates to zero over ω (eq. 6.4).

---

## 3. Sourced goldens (exact numbers, all from the paper)

### 3.1 Jupiter-Callisto system parameters (used throughout as the running example)

| Quantity | Value | Source |
|---|---|---|
| Mass parameter µ (Jupiter-Callisto) | **5.667 × 10⁻⁵** | p.9 text, Fig.5.1 caption |
| Jacobi constant of L₂ for Jupiter-Callisto | **C₂ = 3.00618** | p.4 text |
| Reference for case 2 (particle exterior, C_J close to C₂) | **C_J = 3** (example) | throughout |

### 3.2 1:2 mean-motion resonance hyperbolic fixed point (Fig. 5.2)

    pres = (ω_res, a_res) = (0, (2)^{2/3}) = (0, 1.587)

Eigenvalue of linearized map at pres (eq. 5.2):

    λ = (1/2)[2 + µβγ + √(µβγ(µβγ+4))]

where β = -∂f/∂ω|_{ω=ω_res} and γ = 6π a_res^{5/2}.

### 3.3 Phase space example values (Fig. 5.1, Fig. 5.3)

| Quantity | Value | Source |
|---|---|---|
| Map parameters (Fig. 5.1) | µ = 5.667×10⁻⁵, C_J = 3, ā = 1.35 | Fig.5.1 caption |
| Example trajectory initial a₀ | **1.54** (in units of moon orbital radius) | p.12 text, Fig.5.3 |
| r₂,min of example trajectory (at closest approach to Callisto) | **0.0341** | p.13 text |
| Hill's radius r_h = (µ/3)^{1/3} | **0.0266** | p.13 text |
| Note: particle never enters Hill sphere during 25-orbit migration | r₂,min > r_h by 28% | p.13 |
| ω_max (half-angular-width of max-kick zone, ā=1.35) | **0.01π** | p.13 text |
| Width of phase-space domain at iterate n (conservative lower bound) | α(n) = (0.0628)^n | p.13 |

### 3.4 Reachable-orbit example (Fig. 6.1)

Starting at a₀ = 1.54, map parameters µ = 5.667×10⁻⁵, ā = 1.35:
- C_J = 3.00: reaches [~1.1, ~1.8] after ~25 orbits
- C_J = 2.99: wider reachable zone; faster migration (larger kicks per Fig.6.1b)

C_J = 2.99 < C₂ = 3.00618: **capture regime** — particle can enter Moon realm (ballistic
capture possible).

### 3.5 Rotational Invariant Circle (RIC) example (Fig. 6.3d)

    (µ, C_J, ā) = (5×10⁻⁶, 3, 2.5)  →  RIC found at lower edge of 1:4 resonance island

RIC is an upper transport barrier; no orbit in the connected chaotic zone below it can
reach a above it. RICs exist for µ below a critical value (numerically found ≤ 5×10⁻⁶
for this C_J and ā).

### 3.6 Ballistic capture exit (Fig. 6.2c)

Numerically computed exit region location on Σ_e: C_J = 3.005, µ = 5.667×10⁻⁵.
Exit appears near a ≈ 1.1–1.3, ω/π ≈ 0.2–0.8 (approximately, from ASCII figure).
(Figures-only; exact boundary not tabulated in paper.)

---

## 4. Phase-space structure — the key physical picture

**Connected chaotic zone:** Between resonance island chains, the phase space has a single
connected chaotic zone where lobe dynamics (intersecting stable/unstable manifolds of
resonance hyperbolic points) provide "lanes of fast migration" across many semimajor axis
values with few periapsis passages.

**Resonances as hyperbolic/elliptic pairs:** s:r resonance (particle orbits s times for
every r moon orbits) at a_res = (r/s)^{2/3}. Each resonance has 2(r-s) hyperbolic and
elliptic points alternating. 1:2 resonance at a = (2)^{2/3} ≈ 1.587 has 2 points (one
hyperbolic at ω=0, one elliptic at ω=π — from standard reasoning).

**Capture/escape mechanics (§6):**
- At C_J < C₂: L₂ bottleneck opens; Conley-McGehee tubes (the stable manifold of the
  L₂ Lyapunov orbit) create an "exit" region on Σ_e in (ω, a) space; reaching this
  exit = ballistic capture by the moon
- At K → 0: orbit becomes unbound (escape to infinity)
- RICs provide upper bound on reachable a for small µ

**Time reversal:** if a state at a_n can be reached from a₀ in n orbits, then a state at
a₀ can be reached from a_n in n orbits (mirror trajectory). Implies bidirectional design.

---

## 5. Reuse verdict vs. open threads

### #500 (Keplerian-map genome — SPECULATIVE) — **KEYSTONE: HIGH**

This is THE paper for #500. The Keplerian map (eq. 4.2) directly provides the genome
mechanism for a gravity-assist moon-tour cycler search:

**Genome encoding:** A sequence of n periapsis passages, each described by
    (ω_i, K_i), i = 0, 1, ..., n

parameterizes a multi-resonance gravity-assist trajectory under one moon's influence.
The map propagates the state forward: given (ω_0, K_0) and the moon's µ, each subsequent
state follows without any integration of ODEs. This makes the genome extremely cheap to
evaluate.

**Key parameters for the genome:**
- **C_J** (Jacobi constant): fixed per trajectory; determines which e is associated with
  each a via the Tisserand relation
- **µ** (moon mass parameter): determines the kick magnitude scale (ΔK ~ µ f)
- **ā** (reference K̄ for the kick table): pick ā near the center of the expected a range
- **Kick table f(ω)**: pre-computed once per (C_J, µ, ā); interpolated at runtime

**Capture condition for tour design:**
- Target the exit region on Σ_e (ballistic capture at C_J < C₂)
- The map tells you in O(n) evaluations whether a given (ω_0, K_0) sequence reaches the
  exit in ≤ n passes, without any ODE integration

**Control / search strategy:**
- To decrease a: hold ω ≈ ω_max (just ahead of moon, f < 0)
- To increase a: hold ω ≈ -ω_max (just behind moon, f > 0)
- For a tour: use one moon's map to decrease a to a target value, then switch to the next
  inner moon's map
- Genetic or random search over (ω_0, K_0) sequences is practical: 25 periapsis passages
  cover the migration range shown in Fig.5.3

**Extension to multiple moons (mentioned in §7 as future work):** "extension of the
Keplerian map to include ... (ii) multiple perturbers." This is the natural next step
for the #500 genome: chain moon-by-moon Keplerian maps, each with its own µ and kick
table, switching at the control-transfer point.

**Implementation checklist for #500:**
1. Compute f(ω) lookup table for each moon (µ, C_J, ā) — one quadrature pass per
   table, O(grid_size) evals
2. Implement map eq. (4.2): two-line update per periapsis passage
3. Integrate Tisserand relation to get (a, e) from (K, C_J) at each step
4. Define capture exit: Conley-McGehee tube boundary in (ω, a) — either numerically
   via full CR3BP or analytically approximated by C_J ≤ C₂ and a ≤ some threshold
5. Define escape boundary: K → 0 lobe (eq. 6.3 / Fig.6.3a)
6. Run map forward from initial (ω_0, K_0); test exit/escape after each step

### #318 (joint moon-tour search) — **MODERATE**

The Keplerian map phase space gives the foundation for what joint_cell.py should search
in the multi-resonance inter-moon leg. The connected chaotic zone with resonance island
structure is exactly the space of reachable inter-moon trajectories. The capture exit
region is the "landing condition" for the ballistic-capture building block. A tour
sequence in the Tisserand plane corresponds to a sequence of Keplerian map trajectories,
one per moon leg.

For joint_cell.py / joint_sobol.py: the Keplerian map could replace the direct 3BP
integration for the inter-moon transfer legs, yielding a 2D map search instead of ODE
integration over each leg.

### #465 (multi-rev leveraging) — **HIGH**

The Keplerian map is the analytic quantification of the multi-rev resonant gravity assist
phenomenon. Specifically:
- Fig.5.3 shows how 25 periapsis passages can migrate a from 1.54 to ~1.3 (ΔV ≈ 0)
- The width α(n) = (0.0628)^n of the optimal phase-space domain shrinks rapidly — this
  quantifies how **precisely** the initial ω must be controlled for efficient multi-rev
  leveraging
- For #465: the kick function f and the map update provide the analytical basis for
  predicting how many revolutions are needed to achieve a desired Δa, and what the
  control precision requirement is

### #494 (binary/circumbinary) — **NOT RELEVANT**

Single central body (Jupiter) throughout. No binary or circumbinary content.

---

## 6. Implementation notes

**Kick table:** f(ω) is odd in ω, so only need half the grid. Evaluation via numerical
quadrature of ΔG integral (eq.3.7). The integral is well-behaved for orbit periapsis
outside the Hill sphere (Q > 1, i.e. a > 1 in normalized units).

**Limitation of the approximation (§4):** The map breaks the time-reversal symmetry
of the full PCR3BP, visible as "warping" near ω = 0 in Fig.5.1. For precision design
this matters; for a genome-level screening tool (identifying candidate trajectories for
later full integration) it is acceptable.

**CJ → e coupling:** Must always use the Tisserand / Jacobi relation to get e from
(K, C_J) at each step. Approximation: C_J ≈ 1/a + 2√(a(1-e²)), so
    e² ≈ 1 - [(C_J - 1/a) / 2]² / a
valid when O(µ) corrections are small (satisfied for the Galilean moons µ ≈ 10⁻⁵).

---

## 7. Corpus-index status

**Verdict: `mined`** — full paper read; map equations, phase-space structure, and all
recoverable numerical values extracted. No catalogue rows (purely dynamical methods
paper; no Earth-Mars cycler content). Warrants citation in #500 specification.

---

## 8. Not tabulated / figures-only

- Fig.3.1a (f vs. ω for C_J=3, a=15): exact shape described qualitatively; max kick
  values for specific a visible only in ASCII distortion
- Fig.5.1 (Keplerian map phase portrait, µ=5.667e-5, C_J=3, ā=1.35): resonance island
  locations visible but island widths not tabulated
- Fig.6.2c (numerically computed capture exit on Σ_e, C_J=3.005): location approximate
  from text description; exact boundary requires own numerical computation
- Fig.6.3d (RIC boundary for µ=5e-6): visible at a ≈ 2.47–2.52 but exact a(ω) not
  given in text
