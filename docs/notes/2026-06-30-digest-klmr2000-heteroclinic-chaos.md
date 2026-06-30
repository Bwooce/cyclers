# Digest: Koon–Lo–Marsden–Ross 2000, Chaos 10(2) 427–469

**File:** `koon-lo-marsden-ross-2000-heteroclinic-connections-resonance-transitions-chaos-10-2.pdf`
**Citation:** W. S. Koon, M. W. Lo, J. E. Marsden, S. D. Ross, "Heteroclinic connections between
periodic orbits and resonance transitions in celestial mechanics," *Chaos* **10**(2), 427–469 (2000).
DOI: 10.1063/1.166509.
**Digested:** 2026-06-30.

---

## 1. Model and Method

Planar circular restricted three-body problem (PCR3BP). Two primary masses m_S = 1−μ and
m_J = μ rotating in a common frame; massless particle (comet / spacecraft). Jacobi constant
C = −2E (sign convention: C = −(ẋ²+ẏ²) + 2Ω(x,y)).

Energy-level labelling: C₁ > C₂ > C₃ > C₄ = C₅. Main interest is Case 3:
C₂ > C > C₃ — both L₁ and L₂ necks open, so transit is possible between interior (Sun), Jupiter,
and exterior regions.

**Linearisation near Li (Sec. II):** Saddle × centre eigenvalues ±λ and ±iν. Eigenvectors u₁, u₂,
w₁, w₂. Manifold tubes separate transit (inside tube) from nontransit (outside). McGehee's
four-hemisphere classification on the bounding sphere. Section convention: U₁ = {y=0, x<0}
interior; U₂ = {x=1−μ, y<0} and U₃ = {x=1−μ, y>0} Jupiter; U₄ = {y=0, x<−1} exterior.

**Construction of heteroclinic connection (Sec. III E):** Take equal-energy L₁ and L₂ Lyapunov
orbits (same Jacobi constant). Propagate L₁ unstable manifold tube W^{u,J}_{L1} and L₂ stable
manifold tube W^{s,J}_{L2} into the Jupiter region J and intersect with Poincaré section x = 1−μ.
Record q-th cut Γ^{u,J,q}_{L1} and p-th cut Γ^{s,J,p}_{L2} as closed curves in (y, ẏ). A point in
Γ^{u,J,q}_{L1} ∩ Γ^{s,J,p}_{L2} is a (q,p)-heteroclinic point; the sum q+p must be even. The minimum
cut order for transversal intersection depends on (μ, C).

**Itinerary construction (Sec. IV E):** Track pre-images and images of manifold tube caps through
successive Poincaré sections to isolate phase-space regions with prescribed finite itineraries. The
Poincaré map P satisfies generalised Conley–Moser conditions (strips map to strips with expansion),
guaranteeing existence of all admissible bi-infinite itineraries near the chain.

---

## 2. Sourced Goldens (numbers in text, not read from figures)

### System parameters used throughout
| Quantity | Value | Source |
|---|---|---|
| μ (Sun–Jupiter) | 0.0009537 | Sec. III D, p. 441 |
| C (Oterma encounters) | 3.03 | Sec. III F, p. 443 |
| C (itinerary example) | 3.038 | Sec. IV E p. 450, "just below C₂" |
| C (homoclinic demo) | 3.037 | Sec. III E, p. 442 |
| μ (homoclinic 1,3 example) | 0.1 | Fig. 21 caption p. 442 |
| ΔC for that example | C₁−C = 0.0743 | Fig. 21 caption |

### Linearisation constants at L₁ and L₂ for μ = 0.0009537 (Sec. II B, p. 432)
| Point | a | b |
|---|---|---|
| L₁ | 9.892 | 3.446 |
| L₂ | 8.246 | 2.623 |

### Heteroclinic intersection (Sec. III E, p. 442)
- μ = 0.0009537, C = 3.037.
- First intersection of Γ^{u,J,2}_{L1} and Γ^{s,J,2}_{L2} on the x = 1−μ section: **minimum (q,p)=(2,2)**.
- Intersection near **y ≈ 0.042** (black dots in Fig. 22b; figure-derived — caption states this
  approximate y-value; flag as figure-derived, not tabulated IC).
- Number of revolutions around Jupiter = (q+p−1)/2 = 1.5 (i.e., 1½ revolutions for a (2,2)-connection).

### Resonance connections at C = 3.03 (Sec. V B/C)
- L₁ manifold first cuts meet in **3:2 resonance** interior (a^{-3/2} ≈ 3/2, i.e. a ≈ (2/3)^{2/3}).
- L₂ manifold first cuts meet at **2:3 resonance** (exterior) and at **1:2 resonance**.
- At L = 1.26, ḡ is an angular variable mod 2π; the two 2:3 intersections near L=1.26 are identified.

**No explicit IC table in the paper.** All specific states are given graphically (Figs. 22, 23, 36, 37,
39–41) and are figure-derived.

---

## 3. Algorithm Detail (most reusable for #314 / #411 / #496)

The method (paraphrased from Sec. III E) for finding a heteroclinic connection:

1. Fix μ and C (same for both orbits — mandatory).
2. Correct each Lyapunov orbit with Newton shooting; extract Floquet eigenvectors (STM monodromy).
3. Seed manifold tubes with ε-perturbations along Floquet directions at evenly spaced phases τ.
4. Propagate unstable branch of L₁ forward, stable branch of L₂ backward into the Jupiter region.
5. Record q-th and p-th intersections with the Poincaré section x = 1−μ as closed curves Γ in (y, ẏ).
6. Find intersection points of these two closed curves — these are the heteroclinic points.
7. Integrate forward+backward from each intersection point to trace the full connection.

This is exactly the procedure implemented in `genome/heteroclinic_cycle.py` (#314) and
(with cross-system modifications) in `genome/cross_system_cycle.py` (#405/#411).

**Key constraint:** Both orbits must share EXACTLY the same Jacobi constant. Small C-mismatch
places tubes on different energy surfaces and prevents natural intersection.

---

## 4. Reuse Assessment

### #314 `heteroclinic_cycle.py` (planar CR3BP Lyapunov connection)
**Direct golden reference.** KLMR 2000 is the canonical paper underlying #314's implementation.
The Sun-Jupiter system at μ=0.0009537, C=3.037, and the (2,2)-heteroclinic intersection are
the natural validation targets. The (2,2) crossing occurring in the Jupiter region x=1−μ section
at y≈0.042 (figure-derived, not tabulated) can be used as a positive-control check. The paper
itself is the reference for Wilczak–Zgliczyński (W-Z) validation logic in heteroclinic_cycle.py.

### #411 cross-system Newton stall at |R|=0.59 rad
**Indirect but informative.** KLMR 2000 is an INTRA-system method (both orbits in one CR3BP
frame, patched on x=1−μ in the rotating frame). The cross-system stall in #411 is a different
problem (two separate frames, inertial patch section, 3-free-variable Newton over τ_u/τ_s/θ).
KLMR 2000 does NOT provide a cross-system correction method — no analog of the inter-frame
bridge exists in this paper. However, it clarifies the cut-counting strategy: if our manifolds are
only using the first-crossing detection, we may need to allow q>1 and p>1 cuts. The stall at
|R|=0.59 rad is more likely a basin/branch issue than a cut-count issue, but multi-cut exploration
(analogous to the (2,2) solution here) is one remediation to try in #496.

### #496 two-phase corrector port
**Provides the Poincaré-map framework.** The Braik-Ross 2025 approach (multi-crossing Poincaré
map, then closest-pair as initial guess for differential correction) is the operational version of what
KLMR 2000 describes at the level of existence theory. KLMR 2000 confirms: for multi-revolution
connections, use higher-order cuts (q,p > 1) and that the nearest-approach geometry is in the
correct region. Both the q-cut strategy and the "first intersection of closed Γ curves" language are
directly reusable as the algorithmic substrate of #496.

### #494 (k₁,k₂) cycler construction (tube cross-section seeding)
**Indirect.** KLMR 2000 does not give a scaling law for tube cross-section size vs ΔC. The tube
cross-section on the Poincaré section is the "cap" d⁺ bounded by the asymptotic circle a⁺; the
paper shows this area goes to zero with energy (Sec. II C, "size of the ellipse goes to zero with E")
but does not quantify the scaling rate. See Ross-BozorgMagham 2018 digest for the explicit linear law.

---

## 5. Corpus Status Recommendation
**`mined`** — already serves as the underlying reference for #314/heteroclinic_cycle.py and
cited as ref [6] in Braik-Ross 2025. No new catalogue rows are mined from this paper. Deserves
a CORPUS_INDEX line with status `mined-by-#314`.
