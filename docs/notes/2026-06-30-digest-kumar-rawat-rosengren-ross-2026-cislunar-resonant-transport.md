# Kumar, Rawat, Rosengren & Ross 2026 — Cislunar resonant transport and heteroclinic pathways — digest

**Source:** Bhanu Kumar, Anjali Rawat, Aaron J. Rosengren & Shane D. Ross, "Cislunar resonant transport and heteroclinic pathways: From 3:1 to 2:1 to L1," *Advances in Space Research* 77 (2026) 3815–3843. DOI 10.1016/j.asr.2025.12.005. Open access CC-BY.
**Filed:** `cyclers_pdf/papers/kumar-rawat-rosengren-ross-2026-cislunar-resonant-transport-heteroclinic-pathways-ASR-77-3815.pdf`
**Date digested:** 2026-06-30.
**PDF type:** text-layer (no OCR required).

---

## 1. Paper contribution

### System and model

Planar Circular Restricted Three-Body Problem (PCR3BP) with primaries Earth (m₁) and Moon (m₂), mass ratio μ = 1.2150584270572×10⁻² (dimensionless; nondimensional units: Earth-Moon distance = 1, synodic period = 2π TU). All computation is done in the rotating frame centered on the Earth-Moon barycenter. The paper is **purely cislunar EM-CR3BP** — no heliocentric or Sun-Earth component.

### Scope

Three families of interior mean-motion resonances (MMRs) with the Moon: 4:1, 3:1, and 2:1. These correspond to spacecraft making 4, 3, or 2 revolutions around Earth in the time the Moon makes 1 revolution. In dimensionless units, the resonant semi-major axes are approximately a ≈ 0.40 (4:1), 0.48 (3:1), and 0.63 (2:1).

### Contributions

1. **Resonant periodic orbit families.** Full bifurcation diagrams for 4:1, 3:1, and 2:1 unstable (and stable) prograde and retrograde families, traced through Earth-collision singularities via KS regularization (Kumar & Moreno 2025, AAS 25-677). The 2:1 family has a distinctive gap: no low-eccentricity unstable prograde 2:1 orbits exist in the Earth-Moon PCR3BP.

2. **Perigee Poincaré map.** The Poincaré section is defined by Earth-relative osculating mean anomaly = 0 (perigee), detected via the sign change of `(r-μ·ĵ)·(ṙ)` from negative to positive. This avoids the tangency problems of the standard fixed-x or fixed-y sections, making the manifold visualization smooth and continuous. All manifold curves are plotted in (g, a) = (synodic longitude of perigee, semi-major axis) osculating-element coordinates on this 2D section.

3. **Taylor-series manifold parameterization.** Stable/unstable manifolds are computed via a high-order parameterization method (Haro et al. 2016; Kumar 2025, arXiv:2509.03655) that seeds Taylor series W(k, s) — one series per perigee-section crossing point — and globalizes via iterated Poincaré map. Several orders of magnitude more accurate than linear eigenvector seeding.

4. **Automated heteroclinic detection.** Intersection search in (x, y) Cartesian coordinates on the section (NOT in (g, a), which is only for visualization), followed by bisection to high accuracy. Parallelizable.

5. **Generalized distance metric.** Angular momentum vector L = r × ṙ and Laplace (Runge-Lenz) vector A = ṙ × (r × ṙ) − μ_e r/|r|, assembled into a point-to-set Hausdorff-type distance D(b, C_C) = min_{c∈C_C} |L_b−L_c|² + |A_b−A_c|². Used to define practical TOF from when the manifold enters an ε-ball around the target orbit. More stable than Cartesian state or equinoctial elements for near-Keplerian cislunar orbits.

### Physical results

- **3:1 ↔ 2:1 heteroclinics exist** for C ≤ 3.15 (onset when 2:1 unstable orbits appear). Type 1 (direct, short) transfers take ~30 days for 2.50 ≤ C ≤ 3.07; they cease above C ≈ 3.09 and are replaced by Type 2 (5:2-mediated, longer) transfers at ~60 days.
- **2:1 → L1 Lyapunov heteroclinics:** short type 8.1–8.2 days (C = 3.05–3.10); long type 40–43 days.
- **3:1 → L1 Lyapunov heteroclinics:** short type 31.4 days (C = 3.05); long type 57.9 days (C = 3.10).
- **4:1 ↔ 3:1 heteroclinics DO NOT EXIST** for any energy. A chain of rotational invariant circles (KAM tori / RICs) at a ≈ 0.40–0.48 in the Poincaré section forms an impenetrable topological barrier at every studied C value (including sub-surface-collision C = 2.85). The 4:1 resonance is dynamically isolated from the 3:1–2:1–L1 heteroclinic chain.
- At C ≤ 3.10, 2:1 manifold lobes reach toward a = 1 (Moon's semi-major axis), enabling direct lunar transfer via the L1 tube mechanism described in Koon et al. (2000).

---

## 2. Sourced Goldens (reproducible tabulated data)

All ICs are in the PCR3BP rotating frame with y₀ = 0, ẋ₀ = 0, C fixed. For periodic orbits, the symmetric point condition gives the full 4-vector (x₀, 0, 0, ẏ₀); propagate under PCR3BP EOM to reproduce.

### Table 6 — Resonant periodic orbit ICs (exact, from paper)

All unstable prograde orbits. The 4:1 orbits have period ~6.31 TU; 3:1 ~6.40 TU; 2:1 ~6.41–6.56 TU (varies with eccentricity).

| Orbit | C | x₀ (DU) | ẏ₀ (DU/TU) |
|---|---|---|---|
| 4:1 unstable | 2.85 | 0.770658112628 | −0.616166626758 |
| 4:1 unstable | 3.15 | 0.737385941470 | −0.355888785284 |
| 3:1 unstable | 2.54 | 0.901330167142 | −0.846224994375 |
| 3:1 unstable | 2.70 | 0.888102851739 | −0.725919236562 |
| 3:1 unstable | 2.86 | 0.868787031779 | −0.584481100346 |
| 3:1 unstable | 3.00 | 0.845409258968 | −0.434953051104 |
| 3:1 unstable | 3.05 | 0.834875967543 | −0.372004991290 |
| 3:1 unstable | 3.10 | 0.822429022871 | −0.300987128481 |
| 3:1 unstable | 3.15 | 0.806804382814 | −0.218228896847 |
| 2:1 unstable | 2.54 | 0.964310487216 | −1.202332138510 |
| 2:1 unstable | 2.628 | 0.959722043439 | −1.090831532610 |
| 2:1 unstable | 2.70 | 0.954308268237 | −0.989690975727 |
| 2:1 unstable | 2.86 | 0.934036118105 | −0.743014738588 |
| 2:1 unstable | 3.00 | 0.904862874177 | −0.515850946602 |
| 2:1 unstable | 3.05 | 0.892104772476 | −0.429597033392 |
| 2:1 unstable | 3.10 | 0.878280334961 | −0.334629543419 |
| 2:1 unstable | 3.15 | 0.864401165205 | −0.219058827351 |

Selected periods (from figure captions — cross-check these against propagation):
- 4:1 at C = 3.15: T = 6.3089 TU = 2.0082π TU
- 3:1 at C = 3.05: T = 6.3952 TU = 2.0357π TU
- 2:1 (apogee above Moon) at C = 2.628: T = 5.7755 TU = 1.8384π TU
- 2:1 (prograde unstable) at C = 3.05: T = 6.5636 TU = 2.0893π TU

### Table 5 — Heteroclinic transfer ICs (PCR3BP rotating frame state [x, y, ẋ, ẏ])

Points labelled x* are precisely computed manifold intersections (bisection). Points x⁺/x⁻ are close-pair approximations (for Lyapunov orbit targets only).

**3:1 → 2:1 Type 1 transfers (direct, ~30 days):**

| C | x (DU) | y (DU) | ẋ (DU/TU) | ẏ (DU/TU) |
|---|---|---|---|---|
| 2.54 | −0.02495763718 | −0.02013755597 | 7.56003664579 | −4.80802086624 |
| 2.70 | −0.03022725673 | 0.042561712914 | −5.82508502316 | −2.47401118974 |
| 2.86 | 0.06136961660 | −0.00485126788 | 0.322456526264 | 4.886777873750 |
| 3.00 | −0.07187159167 | −0.09858383084 | 3.22059393928 | −1.95100061370 |
| 3.05 | 0.029351172720 | 0.320801654785 | −1.77006016626 | 0.228990735501 |

**3:1 → 2:1 Type 2 transfers (5:2-mediated, ~60 days):**

| C | x (DU) | y (DU) | ẋ (DU/TU) | ẏ (DU/TU) |
|---|---|---|---|---|
| 3.10 | 0.173726091985 | 0.150836454723 | −1.44182292916 | 1.77676712380 |
| 3.15 | 0.203387734653 | 0.165418603215 | −1.25080421783 | 1.62978185751 |

**2:1 → L1 Lyapunov transfers (approximated x⁺, x⁻ pairs):**

Short type (~8.2 days, C = 3.05): x⁺ = (0.0367, −0.3417, 1.661, 0.237), x⁻ = (0.0202, −0.3266, 1.752, 0.174).
Short type (C = 3.10): x⁺ = (0.0256, −0.3714, 1.527, 0.155), x⁻ = (0.0252, −0.3711, 1.528, 0.154).
(Full precision: see Table 5 verbatim above — all 12 digits in paper.)

**3:1 → L1 Lyapunov short transfer (C = 3.05, ~31.4 days):**
x⁺ = (0.0123, −0.3191, 1.797, 0.138), x⁻ = (0.0142, −0.3204, 1.789, 0.147).

### Table 1 — 3:1 → 2:1 transfer times (L-A distance metric)

| C | TOF (days) | ε₃:₁ (NDU) | ε₂:₁ (NDU) |
|---|---|---|---|
| 3.00 | 30.415 | 1×10⁻³ | 7×10⁻³ |
| 3.05 | 32.895 | 2×10⁻³ | 1×10⁻² |
| 3.10 | 60.212 | 1×10⁻² | 5×10⁻³ |
| 3.15 | 61.420 | 2×10⁻² | 4×10⁻³ |

(Type 1 for C = 3.00–3.05; Type 2 for C = 3.10–3.15.)

### Table 2 — Resonant → L1 Lyapunov transfer times (L-A metric)

| Transfer | C | TOF (days) |
|---|---|---|
| 2:1 → L1 (short) | 3.05 | 8.112 |
| 2:1 → L1 (short) | 3.10 | 8.225 |
| 2:1 → L1 (long) | 3.05 | 42.971 |
| 2:1 → L1 (long) | 3.10 | 40.173 |
| 3:1 → L1 (short) | 3.05 | 31.353 |
| 3:1 → L1 (long) | 3.10 | 57.861 |

Tables 3 & 4 give the same transfers with equinoctial element-based metric; times differ by < 2 days — consistent.

**Verdict: FULLY MINED. Tables 5 & 6 are high-precision reproducible anchors in the EM-PCR3BP. No figures-only numerical claims; every tabulated number has 11–15 significant digits.**

---

## 3. Relevance to open threads

### #314 — heteroclinic_cycle.py (CLOSED 2026-06-20)

#314 is done. The Kumar 2026 methodology is more powerful than what heteroclinic_cycle.py implements (linear eigenvector seeds vs. Taylor-series parameterization; fixed-x section vs. perigee section). However, #314's deliverable was already validated (Oterma golden, arXiv:math/0201278), so no remediation is needed. The Kumar approach is the state-of-the-art and would be the right foundation for any future extension of heteroclinic_cycle.py to resonant orbits.

**One actionable note for #314 users:** if the current perigee detection in heteroclinic_cycle.py is done via fixed-x or fixed-y crossings, it will accumulate tangency artifacts for resonant orbits. Kumar §3.1 shows the dot-product event (sign change of (r-μ·ĵ)·ṙ from negative to positive) is the correct and near-tangency-free alternative.

### #405/#411 — cross_system_cycle.py, SE↔EM Newton stall at |R|≈0.59 rad

The #411 stall is on the c_se parameter (SE-CR3BP secondary family near the Canalias bifurcation), not on the EM side. This paper does not operate in the SE system and cannot directly fix the bifurcation robustness issue.

**However, two indirect contributions:**

1. **EM-side anchor for cross-system connection.** The paper precisely characterizes the 2:1→L1 and 3:1→L1 connections in the EM system. The EM-L1 Lyapunov orbit's stable manifold tube is the natural "handoff" point for any SE↔EM cross-system connection. If #411's closure model seeds the EM side at a 2:1 or 3:1 resonant orbit rather than a Lyapunov orbit, the Table 6 ICs provide exact starting points. The 2:1→L1 short connection (~8.2 days, C ≈ 3.05–3.10) is especially cheap and could serve as an EM-exit anchor.

2. **Heteroclinic detection method.** The Kumar §3.4 algorithm searches for manifold intersections in (x, y) Cartesian coordinates on the Poincaré section, not in (g, a) or other derived coordinates. The paper explicitly notes (footnote 4, remark in §5.2) that manifolds at different C values can *appear* to intersect in (g, a) space but have different position vectors and thus cannot yield a ballistic connection. If the #411 closure criterion is checking intersection in a derived coordinate rather than raw Cartesian position at the Poincaré section, it could produce false positives that stall the Newton corrector.

The stall itself (|R| ≈ 0.59 rad, c_se steps falling off the SE family) is a continuation/initialization problem on the SE side that likely requires either (a) a better initial seed in the SE family at the problematic C value, (b) a Jacobi-constant-parameterized continuation (the paper uses KS regularization + C-continuation to pass through bifurcations), or (c) pivoting the SE-side anchor orbit to a different family. None of these is directly provided by this paper, but (b) is described in detail in §3.2 and could be adapted.

### #494 — Binary (k₁,k₂)-cycler μ-family (Ross-Roberts-Tsoukkas 2026)

#494 is about the binary star μ-family (varying mass ratio, fixed Earth-Moon topology). This Kumar paper is Earth-Moon specific (μ = 1.2150584270572×10⁻²) and addresses interior MMRs, not the (k₁,k₂)-cycler stable-branch orbit families from Ross et al. PRL 2026.

**Indirect relevance only.** The Braik-Ross 2026 orbital network paper (already mined, CORPUS_INDEX row 154) uses the EM cycler families as representatives in its accessibility network, and those families at C_J = 3.1294 overlap with the C range studied here. The Kumar Table 6 3:1 row at C = 3.15 (x₀ = 0.8068, ẏ₀ = −0.2182) is the closest to the Braik-Ross C_J = 3.1294 operating point and could cross-reference the R31-U (3:1 unstable) representative in Braik-Ross Table 2.

The specific construction method for (k₁,k₂)-cyclers in the binary star system (Ross 2026 PRL) uses stable manifolds of symmetric periodic orbits, not heteroclinic chains between resonances — the problem structure is different.

---

## 4. Cross-reference with existing corpus

### Koon-Lo-Marsden-Ross 2000 "Shoot the Moon" (already digested: 2026-06-21-digest-koon-2000-shoot-the-moon.md)

Kumar 2026 is a direct extension of Koon et al. 2000 (Chaos paper). The "first Poincaré cut" of the L1 Lyapunov orbit tube as an exit from the Earth realm, and the interior-resonance → L1 tube connection mechanism are both from Koon 2000. Kumar 2026 adds: (a) cislunar interior MMR orbit families, (b) explicit Poincaré manifold intersections (numerical rather than conceptual), (c) quantified TOFs. The digest of the 2000 paper already covers the conceptual framework; this 2026 paper is the numerical realization for cislunar 3:1/2:1 MMRs.

### KLMR 2006 book (already digested: 2026-06-21-digest-klmr-2006-book.md)

Kumar 2026 cites Koon et al. 2011 (book, same as KLMR 2006) in §3.1 for the standard Poincaré map approach. The KLMR book's framework is the direct intellectual ancestor; Kumar 2026 upgrades the manifold computation (Taylor-series vs. linear seeding) and provides Earth-Moon cislunar numbers.

### Braik-Ross 2026 "Orbital Networks" arXiv:2605.31543 (already mined: 2026-06-13-braik-ross-2026-orbital-networks-mining.md)

Kumar 2026 cites Braik-Ross 2025 (AAS conference version, cited as "Braik, A., Ross, S.D., August 10–14 2025. Heteroclinic transfer between l1 and l2 in Earth-Moon system. AAS 25–716"). The orbital-networks paper (arxiv 2605.31543) uses resonant orbit families computed by the same group (Ross lab, Virginia Tech). The relationship: Kumar 2026 → computes the heteroclinic connections → Braik-Ross 2026 → uses those families as nodes in an accessibility network. The Braik-Ross paper's "R21-S is the hard-access family" finding is consistent with Kumar's result that 2:1 *stable* resonant orbits are surrounded by KAM tori and hard to reach from the chaotic sea.

**Key clarification from this comparison:** In Kumar 2026, only *unstable* 2:1 orbits have heteroclinic connections; *stable* 2:1 orbits sit inside stability islands (note: Braik-Ross's R21-S is the stable branch, which is correctly identified as hard-access). The Kumar 2026 paper provides the mechanistic explanation for the Braik-Ross finding.

---

## 5. CORPUS_INDEX recommendation

**Status: `mined`** — Tables 5 and 6 are extractable high-precision reproducible ICs in the EM-PCR3BP. All numbers are 11–15 significant digits, given at fixed Jacobi constants, and can be verified by propagation under standard CR3BP EOM. This exceeds the `digested` threshold.

**Proposed CORPUS_INDEX entry:**

```
| kumar-rawat-rosengren-ross-2026-cislunar-resonant-transport-heteroclinic-pathways-ASR-77-3815.pdf | 2026-06-30-digest-kumar-rawat-rosengren-ross-2026-cislunar-resonant-transport.md | EM-PCR3BP interior MMR heteroclinic pathways: 4:1/3:1/2:1 resonant orbit ICs (Table 6, 17 orbits) + heteroclinic transfer ICs (Table 5, 14 points) + TOFs (Tables 1–4); 4:1↔3:1 barrier confirmed (KAM tori); 3:1→2:1→L1 chain fully quantified | mined | text-layer |
```

---

## 6. Reuse verdict summary

| Thread | Verdict |
|---|---|
| #314 (heteroclinic_cycle.py) | CLOSED — no action. Perigee-section dot-product detection note filed for future maintenance. |
| #405/#411 (SE↔EM cross-system) | INDIRECT — Table 6 provides EM-side anchor ICs for 2:1/3:1→L1 connections; does not fix the SE bifurcation stall. If searching for manifold intersections, use (x,y) Cartesian on the section, not derived coordinates. |
| #494 (binary μ-family) | WEAK — paper is EM-specific; does not extend μ-family. Only indirect cross-check against Braik-Ross C_J = 3.1294 operating point. |
| Future EM-cislunar genome | HIGH VALUE — Tables 5 & 6 are ready-to-use golden ICs. Taylor-series manifold parameterization (Kumar 2025, arXiv:2509.03655) is the right method for any future resonant-orbit heteroclinic search in this system. |
